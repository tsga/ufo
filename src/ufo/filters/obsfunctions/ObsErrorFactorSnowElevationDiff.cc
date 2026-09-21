/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorSnowElevationDiff.h"

#include <cmath>

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ioda/distribution/Accumulator.h"
#include "ioda/distribution/Distribution.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorFactorSnowElevationDiff> 
    makerObsErrorFactorSnowElevationDiff_("ObsErrorFactorSnowElevationDiff");

// -----------------------------------------------------------------------------

ObsErrorFactorSnowElevationDiff::ObsErrorFactorSnowElevationDiff(
    const eckit::Configuration &config)
  : invars_() {
  oops::Log::trace() << "ObsErrorFactorSnowElevationDiff constructor" << std::endl;
  oops::Log::debug() << "ObsErrorFactorSnowElevationDiff: config = " << config << std::endl;
  
  // Initialize options
  options_.reset(new ObsErrorFactorSnowElevationDiffParameters());
  options_->deserialize(config);

  // Include observation elevation variable
  const std::string obs_elev_var = options_->obs_elevation_var.value();
  invars_ += Variable(obs_elev_var);

  // Include model elevation variable
  const std::string model_elev_var = options_->model_elevation_var.value();
  invars_ += Variable(model_elev_var);
}

// -----------------------------------------------------------------------------

ObsErrorFactorSnowElevationDiff::~ObsErrorFactorSnowElevationDiff() {
  oops::Log::trace() << "ObsErrorFactorSnowElevationDiff destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorFactorSnowElevationDiff::compute(
    const ObsFilterData & data,
    ioda::ObsDataVector<float> & obserr) const {
  oops::Log::trace() << "ObsErrorFactorSnowElevationDiff compute start" << std::endl;
  
  const float missing = util::missingValue<float>();
  oops::Log::warning() << "using missing value " << missing << std::endl;
  const float h_scale = options_->elevation_scale_m.value();
  const float max_infl_factor = options_->max_inflation_factor.value();

  // If no observations on this processor then nothing to do
  if (data.nlocs() == 0) return;

  // Ensure that only one output variable is expected
  ASSERT(obserr.nvars() == 1);

  // Get dimensions
  size_t nlocs = data.nlocs();

  // Get observation elevation variable
  std::vector<float> ob_elevation(nlocs);
  const std::string obs_elev_var = options_->obs_elevation_var.value();
  data.get(Variable(obs_elev_var), ob_elevation);

  // Get model surface elevation variable
  std::vector<float> model_elevation(nlocs);
  const std::string model_elev_var = options_->model_elevation_var.value();
  data.get(Variable(model_elev_var), model_elevation);

  // Compute inflation factor for each observation
  float dz, inflation_factor, inflation_factor_1;
  int iv = 0;
  int count_exparg_toobig = 0;

  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    // If missing observation or model elevation, set obserror to max/missing 
    if (ob_elevation[iloc] == missing || model_elevation[iloc] == missing) {
      obserr[iv][iloc] =  max_infl_factor;  //missing?
    } else {

      // Compute elevation difference
      dz = std::abs(model_elevation[iloc] - ob_elevation[iloc]);

      // inflation factor: inflated error = obs_err / exp(-1 * dz^2 / (h^2))
      inflation_factor_1 = dz * dz / (h_scale * h_scale);

      if (inflation_factor_1 > 40.0) {  // avoid float overflow in arg of exp()
          obserr[iv][iloc] = max_infl_factor;
          count_exparg_toobig = count_exparg_toobig + 1;   // count_exparg_toobig_accumulator->addTerm(iloc, 1);
      } else {
      
          inflation_factor = std::exp(inflation_factor_1);
          // min max cutoff 
	  inflation_factor = std::max(1.0f, inflation_factor); 
          inflation_factor = std::min(max_infl_factor, inflation_factor);

          obserr[iv][iloc] = inflation_factor;
      }
    }
  }
  //const std::size_t count_exparg_toobig = count_exparg_toobig_accumulator->computeResult();
  if (count_exparg_toobig > 0) {
    oops::Log::warning() << "ObsErrorFactorSnowElevationDiff WARNING: " << count_exparg_toobig
      << " instance(s) of exponential argument too large." << std::endl;
  }
  oops::Log::trace() << "ObsErrorFactorSnowElevationDiff compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorSnowElevationDiff::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
