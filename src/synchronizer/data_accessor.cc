/**
 * @file   data_accessor.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  data accessors constructor functions
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
#include "data_accessor.hh"
#include "fe_engine.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T, bool pack_helper>
void DataAccessor<Element>::packUnpackNodalDataHelper(
    Array<T> & data, CommunicationBuffer & buffer,
    const Array<Element> & elements, const Mesh & mesh) {
  UInt nb_component = data.getNbComponent();
  UInt nb_nodes_per_element = 0;

  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  UInt * conn = nullptr;

  for (auto & el : elements) {
    if (el.type != current_element_type ||
        el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type = el.ghost_type;
      conn = mesh.getConnectivity(el.type, el.ghost_type).storage();
      nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);
    }

    UInt el_offset = el.element * nb_nodes_per_element;
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];
      Vector<T> data_vect(data.storage() + offset_conn * nb_component,
                          nb_component);

      if (pack_helper)
        buffer << data_vect;
      else
        buffer >> data_vect;
    }
  }
}

/* ------------------------------------------------------------------------ */
template <typename T, bool pack_helper>
void DataAccessor<Element>::packUnpackElementalDataHelper(
    ElementTypeMapArray<T> & data_to_pack, CommunicationBuffer & buffer,
    const Array<Element> & element, bool per_quadrature_point_data,
    const FEEngine & fem) {
  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  UInt nb_quad_per_elem = 0;
  UInt nb_component = 0;

  Array<T> * vect = nullptr;

  for (auto & el : element) {
    if (el.type != current_element_type ||
        el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type = el.ghost_type;
      vect = &data_to_pack(el.type, el.ghost_type);

      nb_quad_per_elem =
          per_quadrature_point_data
              ? fem.getNbIntegrationPoints(el.type, el.ghost_type)
              : 1;
      nb_component = vect->getNbComponent();
    }

    Vector<T> data(vect->storage() +
                       el.element * nb_component * nb_quad_per_elem,
                   nb_component * nb_quad_per_elem);
    if (pack_helper)
      buffer << data;
    else
      buffer >> data;
  }
}

/* -------------------------------------------------------------------------- */
template <typename T, bool pack_helper>
void DataAccessor<UInt>::packUnpackDOFDataHelper(Array<T> & data,
                                                 CommunicationBuffer & buffer,
                                                 const Array<UInt> & dofs) {
  T * data_ptr = data.storage();
  for (const auto & dof : dofs) {
    if (pack_helper)
      buffer << data_ptr[dof];
    else
      buffer >> data_ptr[dof];
  }
}

/* -------------------------------------------------------------------------- */
#define DECLARE_HELPERS(T)                                                     \
  template void DataAccessor<Element>::packUnpackNodalDataHelper<T, false>(    \
      Array<T> & data, CommunicationBuffer & buffer,                           \
      const Array<Element> & elements, const Mesh & mesh);                     \
  template void DataAccessor<Element>::packUnpackNodalDataHelper<T, true>(     \
      Array<T> & data, CommunicationBuffer & buffer,                           \
      const Array<Element> & elements, const Mesh & mesh);                     \
  template void                                                                \
  DataAccessor<Element>::packUnpackElementalDataHelper<T, false>(              \
      ElementTypeMapArray<T> & data_to_pack, CommunicationBuffer & buffer,     \
      const Array<Element> & element, bool per_quadrature_point_data,          \
      const FEEngine & fem);                                                   \
  template void DataAccessor<Element>::packUnpackElementalDataHelper<T, true>( \
      ElementTypeMapArray<T> & data_to_pack, CommunicationBuffer & buffer,     \
      const Array<Element> & element, bool per_quadrature_point_data,          \
      const FEEngine & fem);                                                   \
  template void DataAccessor<UInt>::packUnpackDOFDataHelper<T, true>(          \
      Array<T> & data, CommunicationBuffer & buffer,                           \
      const Array<UInt> & dofs);                                               \
  template void DataAccessor<UInt>::packUnpackDOFDataHelper<T, false>(         \
      Array<T> & data, CommunicationBuffer & buffer, const Array<UInt> & dofs)

/* -------------------------------------------------------------------------- */
DECLARE_HELPERS(Real);
DECLARE_HELPERS(UInt);
DECLARE_HELPERS(bool);
/* -------------------------------------------------------------------------- */

} // namespace akantu
