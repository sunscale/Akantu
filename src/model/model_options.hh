/**
 * @file   model_options.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Mon Dec 04 2017
 *
 * @brief A Documented file.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_named_argument.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MODEL_OPTIONS_HH__
#define __AKANTU_MODEL_OPTIONS_HH__

namespace akantu {

namespace {
  DECLARE_NAMED_ARGUMENT(analysis_method);
}

struct ModelOptions {
  explicit ModelOptions(AnalysisMethod analysis_method = _static)
      : analysis_method(analysis_method) {}

  template <typename... pack>
  ModelOptions(use_named_args_t, pack &&... _pack)
      : ModelOptions(OPTIONAL_NAMED_ARG(analysis_method, _static)) {}

  virtual ~ModelOptions() = default;

  AnalysisMethod analysis_method;
};

#ifdef AKANTU_SOLID_MECHANICS
/* -------------------------------------------------------------------------- */
struct SolidMechanicsModelOptions : public ModelOptions {
  explicit SolidMechanicsModelOptions(
      AnalysisMethod analysis_method = _explicit_lumped_mass)
      : ModelOptions(analysis_method) {}

  template <typename... pack>
  SolidMechanicsModelOptions(use_named_args_t, pack &&... _pack)
      : SolidMechanicsModelOptions(
            OPTIONAL_NAMED_ARG(analysis_method, _explicit_lumped_mass)) {}
};
#endif

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_COHESIVE_ELEMENT
namespace {
  DECLARE_NAMED_ARGUMENT(is_extrinsic);
}
/* -------------------------------------------------------------------------- */
struct SolidMechanicsModelCohesiveOptions : public SolidMechanicsModelOptions {
  SolidMechanicsModelCohesiveOptions(
      AnalysisMethod analysis_method = _explicit_lumped_mass,
      bool extrinsic = false)
      : SolidMechanicsModelOptions(analysis_method), extrinsic(extrinsic) {}

  template <typename... pack>
  SolidMechanicsModelCohesiveOptions(use_named_args_t, pack &&... _pack)
      : SolidMechanicsModelCohesiveOptions(
            OPTIONAL_NAMED_ARG(analysis_method, _explicit_lumped_mass),
            OPTIONAL_NAMED_ARG(is_extrinsic, false)) {}

  bool extrinsic{false};
};
#endif

#ifdef AKANTU_HEAT_TRANSFER
/* -------------------------------------------------------------------------- */
struct HeatTransferModelOptions : public ModelOptions {
  explicit HeatTransferModelOptions(
      AnalysisMethod analysis_method = _explicit_lumped_mass)
      : ModelOptions(analysis_method) {}

  template <typename... pack>
  HeatTransferModelOptions(use_named_args_t, pack &&... _pack)
      : HeatTransferModelOptions(
            OPTIONAL_NAMED_ARG(analysis_method, _explicit_lumped_mass)) {}
};
#endif

} // akantu

#endif /* __AKANTU_MODEL_OPTIONS_HH__ */
