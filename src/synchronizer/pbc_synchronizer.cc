/**
 * @file   pbc_synchronizer.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Mon Jun 20 2011
 * @date last modification: Tue Nov 06 2012
 *
 * @brief  implementation of the synchronizer for the pbc conditions
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
#include "pbc_synchronizer.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
PBCSynchronizer::PBCSynchronizer(std::map<UInt,UInt> & pairs,
				 SynchronizerID id,
				 MemoryID memory_id):
  Synchronizer(id,memory_id),pbc_pair(pairs){
}
/* -------------------------------------------------------------------------- */
void PBCSynchronizer::asynchronousSynchronize(DataAccessor & data_accessor,
					      SynchronizationTag tag){

  if (size_buffer.find(tag) == size_buffer.end())
    computeBufferSize(data_accessor, tag);

  buffer.resize(size_buffer[tag]);

  for (std::map<UInt,UInt>::iterator it = pbc_pair.begin();
       it != pbc_pair.end();
       ++it) {
    UInt node_master = it->second;
    UInt node_slave = it->first;

    AKANTU_DEBUG_INFO("packing node " << node_master);
    data_accessor.packData(buffer, node_master, tag);

    AKANTU_DEBUG_INFO("unpacking node " << node_slave);
    data_accessor.unpackData(buffer, node_slave, tag);
  }
}

/* -------------------------------------------------------------------------- */
void PBCSynchronizer::waitEndSynchronize(__attribute__((unused)) DataAccessor & data_accessor,
					 __attribute__((unused)) SynchronizationTag tag){
}

/* -------------------------------------------------------------------------- */
void PBCSynchronizer::computeBufferSize(DataAccessor & data_accessor,
					SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  UInt size = 0;

  for (std::map<UInt,UInt>::iterator it = pbc_pair.begin();
       it != pbc_pair.end();++it) {
    size += data_accessor.getNbDataToPack(tag);
  }

  size_buffer[tag] = size;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
