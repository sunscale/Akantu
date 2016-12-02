/**
 * @file   global_ids_updater_inline_impl.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Oct 02 2015
 * @date last modification: Sun Oct 04 2015
 *
 * @brief  Implementation of the inline functions of GlobalIdsUpdater
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
#include "global_ids_updater.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_GLOBAL_IDS_UPDATER_INLINE_IMPL_CC__
#define __AKANTU_GLOBAL_IDS_UPDATER_INLINE_IMPL_CC__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline UInt GlobalIdsUpdater::getNbData(const Array<Element> & elements,
                                        const SynchronizationTag & tag) const {
  UInt size = 0;
  if (elements.getSize() == 0)
    return size;

  if (tag == _gst_giu_global_conn)
    size += Mesh::getNbNodesPerElementList(elements) * sizeof(UInt);

  return size;
}

/* -------------------------------------------------------------------------- */
inline void GlobalIdsUpdater::packData(CommunicationBuffer & buffer,
                                       const Array<Element> & elements,
                                       const SynchronizationTag & tag) const {
  if (tag == _gst_giu_global_conn)
    packUnpackGlobalConnectivity<true>(buffer, elements);
}

/* -------------------------------------------------------------------------- */
inline void GlobalIdsUpdater::unpackData(CommunicationBuffer & buffer,
                                         const Array<Element> & elements,
                                         const SynchronizationTag & tag) {
  if (tag == _gst_giu_global_conn)
    packUnpackGlobalConnectivity<false>(buffer, elements);
}

/* -------------------------------------------------------------------------- */
template <bool pack_mode>
inline void GlobalIdsUpdater::packUnpackGlobalConnectivity(
    CommunicationBuffer & buffer, const Array<Element> & elements) const {
  AKANTU_DEBUG_IN();

  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;

  Array<UInt>::const_vector_iterator conn_begin;
  UInt nb_nodes_per_elem = 0;
  UInt index;

  Array<UInt> & global_nodes_ids = mesh.getGlobalNodesIds();

  Array<Element>::const_scalar_iterator it = elements.begin();
  Array<Element>::const_scalar_iterator end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;

    if (el.type != current_element_type ||
        el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type = el.ghost_type;

      const Array<UInt> & connectivity =
          mesh.getConnectivity(current_element_type, current_ghost_type);
      nb_nodes_per_elem = connectivity.getNbComponent();
      conn_begin = connectivity.begin(nb_nodes_per_elem);
    }

    /// get element connectivity
    const Vector<UInt> current_conn = conn_begin[el.element];

    /// loop on all connectivity nodes
    for (UInt n = 0; n < nb_nodes_per_elem; ++n) {
      UInt node = current_conn(n);

      if (pack_mode) {
        /// if node is local or master pack its global id, otherwise
        /// dummy data
        if (mesh.isLocalOrMasterNode(node))
          index = global_nodes_ids(node);
        else
          index = UInt(-1);

        buffer << index;
      } else {
        buffer >> index;

        /// update slave nodes' index
        if (index != UInt(-1) && mesh.isSlaveNode(node))
          global_nodes_ids(node) = index;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__

#endif /* __AKANTU_GLOBAL_IDS_UPDATER_INLINE_IMPL_CC__ */
