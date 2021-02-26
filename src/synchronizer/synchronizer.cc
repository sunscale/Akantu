/**
 * @file   synchronizer.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 01 2010
 * @date last modification: Wed Nov 15 2017
 *
 * @brief  implementation of the common part
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
#include "synchronizer.hh"
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#include <functional>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Synchronizer::Synchronizer(const Communicator & comm, const ID & id,
                           MemoryID memory_id)
    : Memory(id, memory_id), communicator(comm) {
  int max_tag = comm.getMaxTag();

  this->hash_id = std::hash<std::string>()(this->getID());
  if (max_tag != 0) {
    this->hash_id = this->hash_id % max_tag;
  }

  this->nb_proc = communicator.getNbProc();
  this->rank = communicator.whoAmI();
}

} // namespace akantu
