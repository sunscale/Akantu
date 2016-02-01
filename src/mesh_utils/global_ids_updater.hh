/**
 * @file   global_ids_updater.hh
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Oct 02 2015
 *
 * @brief  Class that updates the global ids of new nodes that are
 * inserted in the mesh. The functions in this class must be called
 * after updating the node types
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
#ifndef __AKANTU_GLOBAL_IDS_UPDATER_HH__
#define __AKANTU_GLOBAL_IDS_UPDATER_HH__

/* -------------------------------------------------------------------------- */
#include "data_accessor.hh"

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

class GlobalIdsUpdater : public DataAccessor {
public:
  GlobalIdsUpdater(Mesh & mesh, DistributedSynchronizer * synchronizer) :
    mesh(mesh), synchronizer(synchronizer) {}

  /// function to update the global connectivity of new inserted
  /// nodes. It must be called after updating the node types.
  UInt updateGlobalIDs(UInt old_nb_nodes);

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  inline virtual UInt getNbDataForElements(const Array<Element> & elements,
					   SynchronizationTag tag) const;

  inline virtual void packElementData(CommunicationBuffer & buffer,
				      const Array<Element> & elements,
				      SynchronizationTag tag) const;

  inline virtual void unpackElementData(CommunicationBuffer & buffer,
					const Array<Element> & elements,
					SynchronizationTag tag);

  template<bool pack_mode>
  inline void packUnpackGlobalConnectivity(CommunicationBuffer & buffer,
					   const Array<Element> & elements) const;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
private:
  /// Reference to the mesh to update
  Mesh & mesh;

  /// distributed synchronizer to communicate the connectivity
  DistributedSynchronizer * synchronizer;
};

__END_AKANTU__

#include "global_ids_updater_inline_impl.cc"

#endif /* __AKANTU_GLOBAL_IDS_UPDATER_HH__ */
