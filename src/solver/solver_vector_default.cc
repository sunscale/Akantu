/**
 * @file   solver_vector.cc
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
#include "solver_vector_default.hh"
#include "dof_manager_default.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SolverVectorArray::SolverVectorArray(DOFManagerDefault & dof_manager,
                                     const ID & id)
    : SolverVector(dof_manager, id), dof_manager(dof_manager) {}

/* -------------------------------------------------------------------------- */
SolverVectorArray::SolverVectorArray(const SolverVectorArray & vector,
                                     const ID & id)
    : SolverVector(vector, id), dof_manager(vector.dof_manager) {}

/* -------------------------------------------------------------------------- */
SolverVectorDefault::SolverVectorDefault(DOFManagerDefault & dof_manager,
                                         const ID & id)
    : SolverVectorArray(dof_manager, id), vector(0, 1, id + ":vector") {}

/* -------------------------------------------------------------------------- */
SolverVectorDefault::SolverVectorDefault(const SolverVectorDefault & vector,
                                         const ID & id)
    : SolverVectorArray(vector, id), vector(vector.vector, id + ":vector") {}

/* -------------------------------------------------------------------------- */
void SolverVectorDefault::resize() {
  this->vector.resize(dof_manager.getLocalSystemSize(), 0.);
}

/* -------------------------------------------------------------------------- */
void SolverVectorDefault::clear() { this->vector.clear(); }

/* -------------------------------------------------------------------------- */
Array<Real> & SolverVectorDefault::getGlobalVector() { return vector; }

/* -------------------------------------------------------------------------- */
void SolverVectorDefault::setGlobalVector(const Array<Real> & solution) {
  vector.copy(solution);
}

} // namespace akantu
