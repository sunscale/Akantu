/**
 * @file   element_info_per_processor_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Mar 16 2016
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Helper classes to create the distributed synchronizer and distribute
 * a mesh
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_group.hh"
#include "element_info_per_processor.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_INFO_PER_PROCESSOR_TMPL_HH__
#define __AKANTU_ELEMENT_INFO_PER_PROCESSOR_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T, typename BufferType>
void ElementInfoPerProc::fillMeshDataTemplated(BufferType & buffer,
                                               const std::string & tag_name,
                                               UInt nb_component) {

  AKANTU_DEBUG_ASSERT(this->mesh.getNbElement(this->type) == nb_local_element,
                      "Did not got enought informations for the tag "
                          << tag_name << " and the element type " << this->type
                          << ":"
                          << "_not_ghost."
                          << " Got " << nb_local_element << " values, expected "
                          << mesh.getNbElement(this->type));

  mesh.getElementalData<T>(tag_name);
  Array<T> & data = mesh.getElementalDataArrayAlloc<T>(
      tag_name, this->type, _not_ghost, nb_component);

  data.resize(nb_local_element);
  /// unpacking the data, element by element
  for (UInt i(0); i < nb_local_element; ++i) {
    for (UInt j(0); j < nb_component; ++j) {
      buffer >> data(i, j);
    }
  }

  AKANTU_DEBUG_ASSERT(mesh.getNbElement(this->type, _ghost) == nb_ghost_element,
                      "Did not got enought informations for the tag "
                          << tag_name << " and the element type " << this->type
                          << ":"
                          << "_ghost."
                          << " Got " << nb_ghost_element << " values, expected "
                          << mesh.getNbElement(this->type, _ghost));

  Array<T> & data_ghost = mesh.getElementalDataArrayAlloc<T>(
      tag_name, this->type, _ghost, nb_component);
  data_ghost.resize(nb_ghost_element);

  /// unpacking the ghost data, element by element
  for (UInt j(0); j < nb_ghost_element; ++j) {
    for (UInt k(0); k < nb_component; ++k) {
      buffer >> data_ghost(j, k);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename BufferType>
void ElementInfoPerProc::fillMeshData(BufferType & buffer,
                                      const std::string & tag_name,
                                      const MeshDataTypeCode & type_code,
                                      UInt nb_component) {
#define AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA(r, extra_param, elem)          \
  case MeshDataTypeCode::BOOST_PP_TUPLE_ELEM(2, 0, elem): {                    \
    fillMeshDataTemplated<BOOST_PP_TUPLE_ELEM(2, 1, elem)>(buffer, tag_name,   \
                                                           nb_component);      \
    break;                                                                     \
  }

  switch (type_code) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA, ,
                          AKANTU_MESH_DATA_TYPES)
  default:
    AKANTU_ERROR("Could not determine the type of tag" << tag_name << "!");
    break;
  }
#undef AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA
}

/* -------------------------------------------------------------------------- */
template <class CommunicationBuffer>
void ElementInfoPerProc::fillElementGroupsFromBuffer(
    CommunicationBuffer & buffer) {
  AKANTU_DEBUG_IN();

  Element el;
  el.type = type;

  for (auto ghost_type : ghost_types) {
    UInt nb_element = mesh.getNbElement(type, ghost_type);
    el.ghost_type = ghost_type;

    for (UInt e = 0; e < nb_element; ++e) {
      el.element = e;

      std::vector<std::string> element_to_group;
      buffer >> element_to_group;

      AKANTU_DEBUG_ASSERT(e < mesh.getNbElement(type, ghost_type),
                          "The mesh does not have the element " << e);

      for (auto && element : element_to_group) {
        mesh.getElementGroup(element).add(el, false, false);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* __AKANTU_ELEMENT_INFO_PER_PROCESSOR_TMPL_HH__ */
