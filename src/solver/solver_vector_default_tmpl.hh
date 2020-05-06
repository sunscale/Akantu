/**
 * @file   solver_vector_default_tmpl.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Jan 01 2019
 *
 * @brief A Documented file.
 *
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
inline SolverVectorArray::SolverVectorArray(DOFManagerDefault & dof_manager,
                                            const ID & id)
    : SolverVector(dof_manager, id) {}

/* -------------------------------------------------------------------------- */
inline SolverVectorArray::SolverVectorArray(const SolverVectorArray & vector,
                                            const ID & id)
    : SolverVector(vector, id) {}

/* -------------------------------------------------------------------------- */
template <class Array_>
SolverVector & SolverVectorArrayTmpl<Array_>::
operator+(const SolverVector & y) {
  const auto & y_ = aka::as_type<SolverVectorArray>(y);
  this->vector += y_.getVector();

  ++this->release_;
  return *this;
}

/* -------------------------------------------------------------------------- */
template <class Array_>
SolverVector & SolverVectorArrayTmpl<Array_>::
operator=(const SolverVector & y) {
  const auto & y_ = aka::as_type<SolverVectorArray>(y);
  this->vector.copy(y_.getVector());

  this->release_ = y.release();
  return *this;
}

/* -------------------------------------------------------------------------- */
template <class Array_> inline Int SolverVectorArrayTmpl<Array_>::size() const {
  return this->dof_manager.getSystemSize();
}

/* -------------------------------------------------------------------------- */
template <class Array_>
inline Int SolverVectorArrayTmpl<Array_>::localSize() const {
  return dof_manager.getLocalSystemSize();
}

} // namespace akantu

#endif /* __AKANTU_SOLVER_VECTOR_DEFAULT_TMPL_HH__ */
