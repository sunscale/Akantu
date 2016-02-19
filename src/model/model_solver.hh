/**
 * @file   model_solver.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Jul 22 10:53:10 2015
 *
 * @brief  Class regrouping the common solve interface to the different models
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
#include "aka_common.hh"
#include "parsable.hh"
#include "solver_callback.hh"
/* -------------------------------------------------------------------------- */
#include <set>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MODEL_SOLVER_HH__
#define __AKANTU_MODEL_SOLVER_HH__

namespace akantu {
  class Mesh;
  class DOFManager;
  class IntegrationScheme;
}

__BEGIN_AKANTU__

class ModelSolver : public Parsable, public SolverCallback {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ModelSolver(const Mesh & mesh, const ID & id, UInt memory_id);
  virtual ~ModelSolver();

  /// initialize the dof manager based on solver type passed in the input file
  void initDOFManager();
  /// initialize the dof manager based on the used chosen solver type
  void initDOFManager(const ID & solver_type);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve a step using a given pre instantiated time step solver and
  /// nondynamic linear solver
  void solveStep(ID time_step_solver_id = "");

  /// Initialize a time solver that can be used afterwards with its id
  void
  getNewSolver(const ID & solver_id,
               TimeStepSolverType time_step_solver_type,
               NonLinearSolverType non_linear_solver_type = _nls_auto);

  /// set an integration scheme for a given dof and a given solver
  void
  setIntegrationScheme(const ID & solver_id, const ID & dof_id,
                       const IntegrationSchemeType & integration_scheme_type);

  /// set an externally instantiated integration scheme
  void setIntegrationScheme(const ID & solver_id, const ID & dof_id,
                            IntegrationScheme & integration_scheme);

  /* ------------------------------------------------------------------------ */
  /* SolverCallback interface                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// Predictor interface for the callback
  virtual void predictor();

  /// Corrector interface for the callback
  virtual void corrector();

  /// AssembleResidual interface for the callback
  virtual void assembleResidual() = 0;

  /// AssembleJacobian interface for the callback
  virtual void assembleJacobian() = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  DOFManager & getDOFManager() { return *this->dof_manager; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  ID parent_id;
  UInt parent_memory_id;

  /// Underlying mesh
  const Mesh & mesh;

  /// Underlying dof_manager (the brain...)
  DOFManager * dof_manager;

  /// Default time step solver to use
  ID default_solver_id;
};

__END_AKANTU__

#endif /* __AKANTU_MODEL_SOLVER_HH__ */
