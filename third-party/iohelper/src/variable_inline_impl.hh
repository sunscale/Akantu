/**
 * @file   variable_inline_impl.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Jun 04 2013
 * @date last modification: Tue Jun 04 2013
 *
 * @brief  inline implementation of dumper visitor
 *
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "paraview_helper.hh"
#include "dumper_lammps.hh"
#include "dumper_text.hh"

#ifndef IOHELPER_VARIABLE_INLINE_IMPL_HH_
#define IOHELPER_VARIABLE_INLINE_IMPL_HH_

/* -------------------------------------------------------------------------- */
namespace iohelper {

template <class Cont>
inline void Variable<Cont>::accept(Visitor & v) const {
  if (auto * ptr_txt = dynamic_cast<DumperText *>(&v)) {
    ptr_txt->visitVariable(*this);
  }
  /*
    else if (ParaviewHelper * ptr_ph = dynamic_cast<ParaviewHelper*>(&v)){
    ptr_ph->visitVariable(*this);
  } else if (iohelper::DumperLammps<iohelper::bond> * ptr_dlb = dynamic_cast<DumperLammps<bond>*>(&v)) {
    ptr_dlb->visitVariable(*this);
  } else if (DumperLammps<iohelper::atomic> * ptr_dla = dynamic_cast<DumperLammps<atomic>*>(&v)){
    ptr_dla->visitVariable(*this);
    }
  */
}

/* -------------------------------------------------------------------------- */

}

#endif /* IOHELPER_VARIABLE_INLINE_IMPL_HH_ */
