/**
 * @file   cohesive_internal_field_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  implementation of the cohesive internal field
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef AKANTU_COHESIVE_INTERNAL_FIELD_TMPL_HH_
#define AKANTU_COHESIVE_INTERNAL_FIELD_TMPL_HH_

namespace akantu {

template <typename T>
CohesiveInternalField<T>::CohesiveInternalField(const ID & id,
                                                Material & material)
    : InternalField<T>(
          id, material, material.getModel().getFEEngine("CohesiveFEEngine"),
          aka::as_type<MaterialCohesive>(material).getElementFilter()) {
  this->element_kind = _ek_cohesive;
}

template <typename T>
CohesiveInternalField<T>::~CohesiveInternalField() = default;

template <typename T>
void CohesiveInternalField<T>::initialize(UInt nb_component) {
  this->internalInitialize(nb_component);
}

/* -------------------------------------------------------------------------- */
template <typename T>
FacetInternalField<T>::FacetInternalField(const ID & id, Material & material)
    : InternalField<T>(
          id, material, material.getModel().getFEEngine("FacetsFEEngine"),
          aka::as_type<MaterialCohesive>(material).getFacetFilter()) {
  this->spatial_dimension -= 1;
  this->element_kind = _ek_regular;
}

template <typename T> FacetInternalField<T>::~FacetInternalField() = default;

template <typename T>
void FacetInternalField<T>::initialize(UInt nb_component) {
  this->internalInitialize(nb_component);
}

/* -------------------------------------------------------------------------- */
template <>
inline void
ParameterTyped<RandomInternalField<Real, FacetInternalField>>::setAuto(
    const ParserParameter & in_param) {
  Parameter::setAuto(in_param);
  RandomParameter<Real> r = in_param;
  param.setRandomDistribution(r);
}

/* -------------------------------------------------------------------------- */
template <>
inline void
ParameterTyped<RandomInternalField<Real, CohesiveInternalField>>::setAuto(
    const ParserParameter & in_param) {
  Parameter::setAuto(in_param);
  RandomParameter<Real> r = in_param;
  param.setRandomDistribution(r);
}

} // namespace akantu

#endif /* AKANTU_COHESIVE_INTERNAL_FIELD_TMPL_HH_ */
