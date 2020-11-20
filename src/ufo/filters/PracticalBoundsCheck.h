/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PRACTICAL_DIFFERENCECHECK_H_
#define UFO_PRACTICAL_DIFFERENCECHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// PracticalBoundsCheck filter

class PracticalBoundsCheck : public FilterBase,
                        private util::ObjectCounter<PracticalBoundsCheck> {
 public:
  static const std::string classname() {return "ufo::PracticalBoundsCheck";}

  PracticalBoundsCheck(ioda::ObsSpace &, const eckit::Configuration &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~PracticalBoundsCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::bounds;}
//   const Variable ref_;
//   const Variable val_;
};

}  // namespace ufo

#endif  // UFO_PRACTICAL_DIFFERENCECHECK_H_
