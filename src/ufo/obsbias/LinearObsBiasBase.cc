/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <string>

#include "ufo/obsbias/LinearObsBiasBase.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

LinearObsBiasFactory::LinearObsBiasFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in ufo::LinearObsBiasFactory."
                       << std::endl;
    ABORT("Element already registered in ufo::LinearObsBiasFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

LinearObsBiasBase * LinearObsBiasFactory::create(const eckit::Configuration & conf) {
  oops::Log::trace() << "LinearObsBiasBase::create starting" << std::endl;
  if (conf.has("ObsBias.name")) {
    std::string id = "";
    id = conf.getString("ObsBias.name");
    typename std::map<std::string, LinearObsBiasFactory*>::iterator jloc = getMakers().find(id);
    if (jloc == getMakers().end()) {
      oops::Log::error() << id << " does not exist in ufo::LinearObsBiasFactory." << std::endl;
      oops::Log::warning() << "Element " << id
                           << " does not exist in ufo::LinearObsBiasFactory." << std::endl;
      return NULL;
    }
    LinearObsBiasBase * ptr = jloc->second->make(conf);
    oops::Log::trace() << "LinearObsBiasBase::create done" << std::endl;
    return ptr;
  } else {
    return NULL;
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo