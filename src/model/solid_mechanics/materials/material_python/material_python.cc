/**
 * @file   material_python.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Fri Nov 13 2015
 * @date last modification: Fri Nov 13 2015
 *
 * @brief  Material python implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "material_python.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
MaterialPython::MaterialPython(SolidMechanicsModel & model, PyObject * obj,
                               const ID & id)
    : Material(model, id), PythonFunctor(obj) {
  AKANTU_DEBUG_IN();

  this->registerInternals();

  std::vector<std::string> param_names =
      this->callFunctor<std::vector<std::string> >("registerParam");

  for (UInt i = 0; i < param_names.size(); ++i) {
    std::stringstream sstr;
    sstr << "PythonParameter" << i;
    this->registerParam(param_names[i], local_params[param_names[i]], 0.,
                        _pat_parsable, sstr.str());
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialPython::registerInternals() {
  std::vector<std::string> internal_names =
      this->callFunctor<std::vector<std::string> >("registerInternals");

  std::vector<UInt> internal_sizes;

  try {
    internal_sizes =
        this->callFunctor<std::vector<UInt> >("registerInternalSizes");
  } catch (...) {
    internal_sizes.assign(internal_names.size(), 1);
  }

  for (UInt i = 0; i < internal_names.size(); ++i) {
    std::stringstream sstr;
    sstr << "PythonInternal" << i;
    this->internals[internal_names[i]] =
        new InternalField<Real>(internal_names[i], *this);
    AKANTU_DEBUG_INFO("alloc internal " << internal_names[i] << " "
                                        << this->internals[internal_names[i]]);

    this->internals[internal_names[i]]->initialize(internal_sizes[i]);
  }

  // making an internal with the quadrature points coordinates
  this->internals["quad_coordinates"] =
      new InternalField<Real>("quad_coordinates", *this);
  auto && coords = *this->internals["quad_coordinates"];
  coords.initialize(this->getSpatialDimension());
}

/* -------------------------------------------------------------------------- */
void MaterialPython::initMaterial() {
  AKANTU_DEBUG_IN();

  Material::initMaterial();

  auto && coords = *this->internals["quad_coordinates"];
  this->model->getFEEngine().computeIntegrationPointsCoordinates(
      coords, &this->element_filter);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialPython::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  std::map<std::string, Array<Real> *> internal_arrays;
  for (auto & i : this->internals) {
    auto & array = (*i.second)(el_type, ghost_type);
    auto & name = i.first;
    internal_arrays[name] = &array;
  }

  auto params = local_params;
  params["rho"] = this->rho;

  this->callFunctor<void>("computeStress", this->gradu(el_type, ghost_type),
                          this->stress(el_type, ghost_type), internal_arrays,
                          params);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename it_type>
void MaterialPython::computeStress(Matrix<Real> & grad_u, Matrix<Real> & sigma,
                                   std::vector<it_type> & internal_iterators) {
  std::vector<Real> inputs;
  for (auto & i : internal_iterators) {
    inputs.push_back(*i);
  }

  for (UInt i = 0; i < inputs.size(); ++i) {
    *internal_iterators[i] = inputs[i];
  }
}

/* -------------------------------------------------------------------------- */
void MaterialPython::computeTangentModuli(const ElementType & el_type,
                                          Array<Real> & tangent_matrix,
                                          GhostType ghost_type) {
  std::map<std::string, Array<Real> *> internal_arrays;
  for (auto & i : this->internals) {
    auto & array = (*i.second)(el_type, ghost_type);
    auto & name = i.first;
    internal_arrays[name] = &array;
  }

  auto params = local_params;
  params["rho"] = this->rho;

  this->callFunctor<void>("computeTangentModuli",
                          this->gradu(el_type, ghost_type), tangent_matrix,
                          internal_arrays, params);
}

/* -------------------------------------------------------------------------- */
Real MaterialPython::getPushWaveSpeed(const Element &) const {
  auto params = local_params;
  params["rho"] = this->rho;

  return this->callFunctor<Real>("getPushWaveSpeed", params);
}

}  // akantu

__END_AKANTU__
