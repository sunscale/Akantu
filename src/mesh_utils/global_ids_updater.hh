/**
 * @file   global_ids_updater.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Oct 02 2015
 * @date last modification: Fri Dec 08 2017
 *
 * @brief  Class that updates the global ids of new nodes that are
 * inserted in the mesh. The functions in this class must be called
 * after updating the node types
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef __AKANTU_GLOBAL_IDS_UPDATER_HH__
#define __AKANTU_GLOBAL_IDS_UPDATER_HH__

/* -------------------------------------------------------------------------- */
#include "data_accessor.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
class ElementSynchronizer;
} // namespace akantu

namespace akantu {

class GlobalIdsUpdater : public DataAccessor<Element> {
public:
  GlobalIdsUpdater(Mesh & mesh, ElementSynchronizer & synchronizer)
      : mesh(mesh), synchronizer(synchronizer) {}

  /// function to update and synchronize the global connectivity of
  /// new inserted nodes. It must be called after updating the node
  /// types. (It calls in sequence the functions
  /// updateGlobalIDsLocally and synchronizeGlobalIDs)
  UInt updateGlobalIDs(UInt local_nb_new_nodes);

  /// function to update the global connectivity (only locally) of new
  /// inserted nodes. It must be called after updating the node types.
  UInt updateGlobalIDsLocally(UInt local_nb_new_nodes);

  /// function to synchronize the global connectivity of new inserted
  /// nodes among the processors. It must be called after updating the
  /// node types.
  void synchronizeGlobalIDs();

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;
  /* ------------------------------------------------------------------------ */
  template <bool pack_mode>
  inline void
  packUnpackGlobalConnectivity(CommunicationBuffer & buffer,
                               const Array<Element> & elements) const;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
private:
  /// Reference to the mesh to update
  Mesh & mesh;

  /// distributed synchronizer to communicate the connectivity
  ElementSynchronizer & synchronizer;

  /// Tells if a reduction is taking place or not
  bool reduce{false};

  std::unordered_map<UInt, std::vector<std::pair<UInt, NodeFlag>>> nodes_flags;
};

} // namespace akantu

#include "global_ids_updater_inline_impl.hh"

#endif /* __AKANTU_GLOBAL_IDS_UPDATER_HH__ */
