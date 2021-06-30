/**
 * @file   grid_synchronizer.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Synchronizer based on spatial grid
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
#include "aka_common.hh"
#include "element_synchronizer.hh"
#include "synchronizer_registry.hh"

/* -------------------------------------------------------------------------- */

#ifndef AKANTU_GRID_SYNCHRONIZER_HH_
#define AKANTU_GRID_SYNCHRONIZER_HH_

namespace akantu {

class Mesh;
template <class T> class SpatialGrid;

class GridSynchronizer : public ElementSynchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  template <typename E>
  GridSynchronizer(Mesh & mesh, const SpatialGrid<E> & grid,
                   const ID & id = "grid_synchronizer",
                   bool register_to_event_manager = true,
                   EventHandlerPriority event_priority = _ehp_synchronizer);

  template <typename E>
  GridSynchronizer(Mesh & mesh, const SpatialGrid<E> & grid,
                   SynchronizerRegistry & synchronizer_registry,
                   const std::set<SynchronizationTag> & tags_to_register,
                   const ID & id = "grid_synchronizer",
                   bool register_to_event_manager = true,
                   EventHandlerPriority event_priority = _ehp_synchronizer);

  ~GridSynchronizer() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
private:
  /**
   *Create the Grid Synchronizer:
   *Compute intersection and send info to neighbours that will be stored in
   *ghosts elements
   */
  template <typename E>
  void createGridSynchronizer(const SpatialGrid<E> & grid);

protected:
  /// Define the tags that will be used in the send and receive instructions
  enum CommTags {
    SIZE_TAG = 0,
    DATA_TAG = 1,
    ASK_NODES_TAG = 2,
    SEND_NODES_TAG = 3
  };

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

} // namespace akantu

#include "grid_synchronizer_tmpl.hh"

#endif /* AKANTU_GRID_SYNCHRONIZER_HH_ */
