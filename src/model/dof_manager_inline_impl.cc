/**
 * @file   dof_manager_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 12 11:07:01 2015
 *
 * @brief
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_DOF_MANAGER_INLINE_IMPL_CC__
#define __AKANTU_DOF_MANAGER_INLINE_IMPL_CC__

__BEGIN_AKANTU__

inline void DOFManager::extractElementEquationNumber(
    const Array<UInt> & equation_numbers, const Vector<UInt> & connectivity,
    UInt nb_degree_of_freedom, Vector<UInt> & local_equation_number) {
  for (UInt i = 0, ld = 0; i < connectivity.size(); ++i) {
    UInt n = connectivity(i);
    for (UInt d = 0; d < nb_degree_of_freedom; ++d, ++ld) {
      local_equation_number(ld) =
          equation_numbers(n * nb_degree_of_freedom + d);
    }
  }
}

__END_AKANTU__

#endif /* __AKANTU_DOF_MANAGER_INLINE_IMPL_CC__ */
