/**
 * @file   grid_synchronizer.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Oct 03 2011
 * @date last modification: Thu Feb 21 2013
 *
 * @brief  synchronizer based in RegularGrid
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "aka_common.hh"
#include "distributed_synchronizer.hh"

/* -------------------------------------------------------------------------- */


#ifndef __AKANTU_GRID_SYNCHRONIZER_HH__
#define __AKANTU_GRID_SYNCHRONIZER_HH__

__BEGIN_AKANTU__

class Mesh;
template<class T>
class SpatialGrid;

class GridSynchronizer : public DistributedSynchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  GridSynchronizer(Mesh & mesh, const ID & id = "grid_synchronizer", MemoryID memory_id = 0);

public:
  virtual ~GridSynchronizer() { };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  template <class E>
  static GridSynchronizer *
  createGridSynchronizer(Mesh & mesh,
			 const SpatialGrid<E> & grid,
			 SynchronizerID id = "grid_synchronizer",
			 MemoryID memory_id = 0);


protected:
  enum CommTags {
    SIZE_TAG        = 0,
    DATA_TAG        = 1,
    ASK_NODES_TAG   = 2,
    SEND_NODES_TAG  = 3
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


__END_AKANTU__

#endif /* __AKANTU_GRID_SYNCHRONIZER_HH__ */
