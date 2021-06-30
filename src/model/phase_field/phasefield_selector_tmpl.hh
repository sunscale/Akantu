/**
 * @file   phasefield_selector_tmpl.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Sun Mar 1 2020
 * @date last modification: Sun Mar 01 2020
 *
 * @brief  Implementation of the template PhaseFieldSelector
 *
 * @section LICENSE
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
#include "phasefield_selector.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASEFIELD_SELECTOR_TMPL_HH__
#define __AKANTU_PHASEFIELD_SELECTOR_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <>
inline UInt ElementDataPhaseFieldSelector<std::string>::
operator()(const Element & element) {
  try {
    std::string material_name = this->elementData(element);
    return model.getPhaseFieldIndex(material_name);
  } catch (...) {
    return PhaseFieldSelector::operator()(element);
  }
}

/* -------------------------------------------------------------------------- */
template <>
inline UInt ElementDataPhaseFieldSelector<UInt>::
operator()(const Element & element) {
  try {
    return this->elementData(element) - first_index;
  } catch (...) {
    return PhaseFieldSelector::operator()(element);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
MeshDataPhaseFieldSelector<T>::MeshDataPhaseFieldSelector(
    const std::string & name, const PhaseFieldModel & model,
    UInt first_index)
    : ElementDataPhaseFieldSelector<T>(model.getMesh().getData<T>(name), model,
                                     first_index) {}

} // namespace akantu

#endif /* __AKANTU_PHASEFIELD_SELECTOR_TMPL_HH__ */
