/**
 * @file   solver_vector_default_tmpl.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Jan 01 2019
 *
 * @brief A Documented file.
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
#include "dof_manager_default.hh"
#include "solver_vector_default.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLVER_VECTOR_DEFAULT_TMPL_HH__
#define __AKANTU_SOLVER_VECTOR_DEFAULT_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Array>
SolverVectorDefaultWrap<Array>::SolverVectorDefaultWrap(
    DOFManagerDefault & dof_manager, Array & vector)
    : SolverVectorArray(dof_manager, vector.getID()), vector(vector) {}

/* -------------------------------------------------------------------------- */
template <class Array>
SolverVectorDefaultWrap<Array>::SolverVectorDefaultWrap(
    const SolverVectorDefaultWrap & vector, const ID & id)
    : SolverVectorArray(vector, id), vector(vector.vector) {}

/* -------------------------------------------------------------------------- */
template <class Array> void SolverVectorDefaultWrap<Array>::resize() {
  static_assert(not std::is_const<Array>::value, "Cannot resize a const Array");
  this->vector.resize(dof_manager.getLocalSystemSize(), 0.);
}

/* -------------------------------------------------------------------------- */
template <class Array> void SolverVectorDefaultWrap<Array>::clear() {
  static_assert(not std::is_const<Array>::value, "Cannot clear a const Array");
  this->vector.clear();
}

inline Int SolverVectorArray::size() { return dof_manager.getSystemSize(); }
inline Int SolverVectorArray::localSize() { return dof_manager.getLocalSystemSize(); }

} // namespace akantu

#endif /* __AKANTU_SOLVER_VECTOR_DEFAULT_TMPL_HH__ */
