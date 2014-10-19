/**
 * @file   pbc_synchronizer.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Jun 16 2011
 * @date last modification: Tue Nov 06 2012
 *
 * @brief  Dofs Synchronizer for periodic boundary condition
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_PBC_SYNCHRONIZER_HH__
#define __AKANTU_PBC_SYNCHRONIZER_HH__
/* -------------------------------------------------------------------------- */
#include "synchronizer.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
/* -------------------------------------------------------------------------- */


class PBCSynchronizer : public Synchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  PBCSynchronizer(std::map<UInt,UInt> & pairs,
		  SynchronizerID id = "pbc_synch",
		  MemoryID memory_id = 0);
  virtual ~PBCSynchronizer(){};

public:

  /* ------------------------------------------------------------------------ */
  /* Inherited from Synchronizer                                              */
  /* ------------------------------------------------------------------------ */

  /// asynchronous synchronization of ghosts
  void asynchronousSynchronize(DataAccessor & data_accessor,SynchronizationTag tag);

  /// wait end of asynchronous synchronization of ghosts
  void waitEndSynchronize(DataAccessor & data_accessor,SynchronizationTag tag);

private:

  /// compute buffer size for a given tag and data accessor
  void computeBufferSize(DataAccessor & data_accessor, SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// pairs of nodes slave -> master node
  std::map<UInt,UInt> & pbc_pair;

  /// size of packing buffer for each tag
  std::map<SynchronizationTag, UInt> size_buffer;

  /// buffer to pack and unpack the data
  CommunicationBuffer buffer;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "pbc_synchronizer_inline_impl.cc"

// /// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const PBCSynchronizer & _this)
// {
//   _this.printself(stream);
//   return stream;
// }

__END_AKANTU__

#endif /* __AKANTU_PBC_SYNCHRONIZER_HH__ */
