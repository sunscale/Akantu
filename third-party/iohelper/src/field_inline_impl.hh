/**
 * @file   field_inline_impl.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Wed Oct 31 2012
 * @date last modification: Tue Jun 04 2013
 *
 * @brief  inline implementation of dumper visitor
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

/* -------------------------------------------------------------------------- */
#include "paraview_helper.hh"
#include "dumper_lammps.hh"
#include "dumper_text.hh"

#ifndef IOHELPER_FIELD_INLINE_IMPL_HH_
#define IOHELPER_FIELD_INLINE_IMPL_HH_

/* -------------------------------------------------------------------------- */
namespace iohelper {

template <class Cont>
inline void Field<Cont>::accept(Visitor & v){
  if (auto * ptr_ph = dynamic_cast<ParaviewHelper *>(&v)) {
    ptr_ph->visitField(*this);
  } else if (auto * ptr_dlb = dynamic_cast<DumperLammps<bond> *>(&v)) {
    ptr_dlb->visitField(*this);
  } else if (auto * ptr_dla = dynamic_cast<DumperLammps<atomic> *>(&v)) {
    ptr_dla->visitField(*this);
  } else if (auto * ptr_txt = dynamic_cast<DumperText *>(&v)) {
    ptr_txt->visitField(*this);
  }
}

}
/* -------------------------------------------------------------------------- */

#endif /* IOHELPER_FIELD_INLINE_IMPL_HH_ */
