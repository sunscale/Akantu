/**
 * @file   solid_mechanics_model_inline_impl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Tue Dec 05 2017
 *
 * @brief  Implementation of the inline functions of the SolidMechanicsModel
 * class
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_named_argument.hh"
#include "material_selector.hh"
#include "material_selector_tmpl.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLID_MECHANICS_MODEL_INLINE_IMPL_HH_
#define AKANTU_SOLID_MECHANICS_MODEL_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline decltype(auto) SolidMechanicsModel::getMaterials() {
  return make_dereference_adaptor(materials);
}

/* -------------------------------------------------------------------------- */
inline decltype(auto) SolidMechanicsModel::getMaterials() const {
  return make_dereference_adaptor(materials);
}

/* -------------------------------------------------------------------------- */
inline Material & SolidMechanicsModel::getMaterial(UInt mat_index) {
  AKANTU_DEBUG_ASSERT(mat_index < materials.size(),
                      "The model " << id << " has no material no "
                                   << mat_index);
  return *materials.at(mat_index);
}

/* -------------------------------------------------------------------------- */
inline const Material & SolidMechanicsModel::getMaterial(UInt mat_index) const {
  AKANTU_DEBUG_ASSERT(mat_index < materials.size(),
                      "The model " << id << " has no material no "
                                   << mat_index);
  return *materials.at(mat_index);
}

/* -------------------------------------------------------------------------- */
inline Material & SolidMechanicsModel::getMaterial(const std::string & name) {
  std::map<std::string, UInt>::const_iterator it =
      materials_names_to_id.find(name);
  if(it == materials_names_to_id.end()) {
    AKANTU_SILENT_EXCEPTION("The model " << id << " has no material named " << name);
  }

  return *materials[it->second];
}

/* -------------------------------------------------------------------------- */
inline UInt
SolidMechanicsModel::getMaterialIndex(const std::string & name) const {
  auto it = materials_names_to_id.find(name);
  if(it == materials_names_to_id.end()) {
    AKANTU_SILENT_EXCEPTION("The model " << id << " has no material named " << name);
  }

  return it->second;
}

/* -------------------------------------------------------------------------- */
inline const Material &
SolidMechanicsModel::getMaterial(const std::string & name) const {
  auto it = materials_names_to_id.find(name);
  if(it == materials_names_to_id.end()) {
    AKANTU_SILENT_EXCEPTION("The model " << id << " has no material named " << name);
  }
  return *materials[it->second];
}

/* -------------------------------------------------------------------------- */
} // namespace akantu

#endif /* AKANTU_SOLID_MECHANICS_MODEL_INLINE_IMPL_HH_ */
