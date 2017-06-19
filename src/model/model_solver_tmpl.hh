/**
 * @file   model_solver_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Apr 29 17:16:33 2016
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
#include "time_step_solver.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MODEL_SOLVER_TMPL_HH__
#define __AKANTU_MODEL_SOLVER_TMPL_HH__

namespace akantu {

template<typename T>
void ModelSolver::set(const ID & param, const T & value, const ID & solver_id) {
  TimeStepSolver & solver = this->getSolver(solver_id);
  solver.setParam(param, value);
}


inline Parameter ModelSolver::get(const ID & param, const ID & solver_id) {
  TimeStepSolver & solver = this->getSolver(solver_id);
  return solver.getParam(param);
}


} // akantu

#endif /* __AKANTU_MODEL_SOLVER_TMPL_HH__ */
