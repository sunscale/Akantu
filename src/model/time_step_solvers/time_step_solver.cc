/**
 * @file   time_step_solver.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Oct 12 16:56:43 2015
 *
 * @brief  Implementation of common part of TimeStepSolvers
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "time_step_solver.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
TimeStepSolver::TimeStepSolver(DOFManager & dof_manager,
                               const TimeStepSolverType & type, const ID & id,
                               UInt memory_id)
    : Memory(id, memory_id), _dof_manager(dof_manager), type(type) {}

/* -------------------------------------------------------------------------- */
TimeStepSolver::~TimeStepSolver() {}

/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */

__END_AKANTU__
