/**
 * @file   solver_vector_distributed.hh
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
#include "solver_vector_default.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLVER_VECTOR_DISTRIBUTED_HH_
#define AKANTU_SOLVER_VECTOR_DISTRIBUTED_HH_

namespace akantu {

class SolverVectorDistributed : public SolverVectorDefault {
public:
  SolverVectorDistributed(DOFManagerDefault & dof_manager,
                          const ID & id = "solver_vector_mumps");

  SolverVectorDistributed(const SolverVectorDefault & vector,
                          const ID & id = "solver_vector_mumps");

  Array<Real> & getGlobalVector() override;
  void setGlobalVector(const Array<Real> & solution) override;

protected:
  // full vector in case it needs to be centralized on master
  std::unique_ptr<Array<Real>> global_vector;
};

} // namespace akantu

#endif /* AKANTU_SOLVER_VECTOR_DISTRIBUTED_HH_ */
