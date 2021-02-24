/**
 * @file   model_solver.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Class regrouping the common solve interface to the different models
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "integration_scheme.hh"
#include "parsable.hh"
#include "solver_callback.hh"
#include "synchronizer_registry.hh"
/* -------------------------------------------------------------------------- */
#include <set>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MODEL_SOLVER_HH__
#define __AKANTU_MODEL_SOLVER_HH__

namespace akantu {
class Mesh;
class DOFManager;
class TimeStepSolver;
class NonLinearSolver;
struct ModelSolverOptions;
} // namespace akantu

namespace akantu {

class ModelSolver : public Parsable,
                    public SolverCallback,
                    public SynchronizerRegistry {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ModelSolver(Mesh & mesh, const ModelType & type, const ID & id,
              UInt memory_id);
  ~ModelSolver() override;

  /// initialize the dof manager based on solver type passed in the input file
  void initDOFManager();
  /// initialize the dof manager based on the used chosen solver type
  void initDOFManager(const ID & solver_type);

protected:
  /// initialize the dof manager based on the used chosen solver type
  void initDOFManager(const ParserSection & section, const ID & solver_type);

public:
  /// Callback for the model to instantiate the matricees when needed
  virtual void initSolver(TimeStepSolverType /*time_step_solver_type*/,
                          NonLinearSolverType /*non_linear_solver_type*/) {}

  /// get the section in the input file (if it exsits) corresponding to this
  /// model
  std::tuple<ParserSection, bool> getParserSection();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve a step using a given pre instantiated time step solver and
  /// nondynamic linear solver
  virtual void solveStep(const ID & solver_id = "");

  /// solve a step using a given pre instantiated time step solver and
  /// non linear solver with a user defined callback instead of the
  /// model itself /!\ This can mess up everything
  virtual void solveStep(SolverCallback & callback, const ID & solver_id = "");

  /// Initialize a time solver that can be used afterwards with its id
  void getNewSolver(
      const ID & solver_id, TimeStepSolverType time_step_solver_type,
      NonLinearSolverType non_linear_solver_type = NonLinearSolverType::_auto);

  /// set an integration scheme for a given dof and a given solver
  void
  setIntegrationScheme(const ID & solver_id, const ID & dof_id,
                       const IntegrationSchemeType & integration_scheme_type,
                       IntegrationScheme::SolutionType solution_type =
                           IntegrationScheme::_not_defined);

  /// set an externally instantiated integration scheme
  void setIntegrationScheme(const ID & solver_id, const ID & dof_id,
                            IntegrationScheme & integration_scheme,
                            IntegrationScheme::SolutionType solution_type =
                                IntegrationScheme::_not_defined);

  /* ------------------------------------------------------------------------ */
  /* SolverCallback interface                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// Predictor interface for the callback
  void predictor() override;

  /// Corrector interface for the callback
  void corrector() override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// Default time step solver to instantiate for this model
  virtual TimeStepSolverType getDefaultSolverType() const;

  /// Default configurations for a given time step solver
  virtual ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const;

  /// get access to the internal dof manager
  DOFManager & getDOFManager() { return *this->dof_manager; }

  /// get the time step of a given solver
  Real getTimeStep(const ID & solver_id = "") const;
  /// set the time step of a given solver
  virtual void setTimeStep(Real time_step, const ID & solver_id = "");

  /// set the parameter 'param' of the solver 'solver_id'
  // template <typename T>
  // void set(const ID & param, const T & value, const ID & solver_id = "");

  /// get the parameter 'param' of the solver 'solver_id'
  // const Parameter & get(const ID & param, const ID & solver_id = "") const;

  /// answer to the question "does the solver exists ?"
  bool hasSolver(const ID & solver_id) const;

  /// changes the current default solver
  void setDefaultSolver(const ID & solver_id);

  /// is a default solver defined
  bool hasDefaultSolver() const;

  /// is an integration scheme set for a given solver and a given dof
  bool hasIntegrationScheme(const ID & solver_id, const ID & dof_id) const;

  TimeStepSolver & getTimeStepSolver(const ID & solver_id = "");
  NonLinearSolver & getNonLinearSolver(const ID & solver_id = "");

  const TimeStepSolver & getTimeStepSolver(const ID & solver_id = "") const;
  const NonLinearSolver & getNonLinearSolver(const ID & solver_id = "") const;

private:
  TimeStepSolver & getSolver(const ID & solver_id);
  const TimeStepSolver & getSolver(const ID & solver_id) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  ModelType model_type;

private:
  ID parent_id;
  UInt parent_memory_id;

  /// Underlying mesh
  Mesh & mesh;

  /// Underlying dof_manager (the brain...)
  std::unique_ptr<DOFManager> dof_manager;

  /// Default time step solver to use
  ID default_solver_id;
};

struct ModelSolverOptions {
  NonLinearSolverType non_linear_solver_type;
  std::map<ID, IntegrationSchemeType> integration_scheme_type;
  std::map<ID, IntegrationScheme::SolutionType> solution_type;
};

} // namespace akantu

#endif /* __AKANTU_MODEL_SOLVER_HH__ */
