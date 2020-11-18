/**
 * @file   dumper_restart.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Thu Dec 06 2012
 *
 * @brief  implementation for restart dumper
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * IOHelper is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include <iomanip>


#include "dumper_restart.hh"
/* -------------------------------------------------------------------------- */
#if defined(__INTEL_COMPILER)
/// remark #981: operands are evaluated in unspecified order
#pragma warning ( disable : 981 )
#endif //defined(__INTEL_COMPILER)

namespace iohelper {

/* -------------------------------------------------------------------------- */


void DumperRestart::dump(__attribute__((unused)) const std::string & basename){

  /* nodal data */
  std::map<std::string,FieldInterface *>::iterator it = per_node_data.begin();
  std::map<std::string,FieldInterface *>::iterator end = per_node_data.end();

  for (; it != end ; ++it) {
    (*it).second->accept(*this);
  }

  /* element data */
  it = per_element_data.begin();
  end = per_element_data.end();
  for (; it != end ; ++it) {
    (*it).second->accept(*this);
  }

  this->incDumpStep();
}
/* -------------------------------------------------------------------------- */





}


/* -------------------------------------------------------------------------- */



