/**
 * @file   model_solver.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 12 13:31:56 2015
 *
 * @brief  Implementation of ModelSolver
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
#include "model_solver.hh"
#include "dof_manager.hh"
#include "non_linear_solver.hh"
#include "time_step_solver.hh"

#include "mesh.hh"

#if defined(AKANTU_USE_MPI)
#include "mpi_type_wrapper.hh"
#endif

#include "dof_manager_default.hh"

#if defined(AKANTU_USE_PETSC)
#include "dof_manager_petsc.hh"
#endif

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ModelSolver::ModelSolver(Mesh & mesh, const ID & id, UInt memory_id)
    : Parsable(_st_model_solver, id), SolverCallback(), parent_id(id),
      parent_memory_id(memory_id), mesh(mesh), dof_manager(NULL),
      default_solver_id("") {}

/* -------------------------------------------------------------------------- */
ModelSolver::~ModelSolver() { delete this->dof_manager; }

/* -------------------------------------------------------------------------- */
void ModelSolver::initDOFManager() {
  std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
    sub_sect = getStaticParser().getSubSections(_st_model_solver);

  // default without external solver activated at compilation same as mumps that
  // is the historical solver but with only the lumped solver
  ID solver_type = "explicit";

#if defined(AKANTU_USE_MUMPS)
  solver_type = "mumps";
#elif defined(AKANTU_USE_PETSC)
  solver_type = "petsc";
#endif

  const ParserSection * section = NULL;
  Parser::const_section_iterator it;
  for (it = sub_sect.first; it != sub_sect.second && section == NULL; ++it) {
    if (it->getName() == this->parent_id) {
      section = &(*it);
      solver_type = section->getOption(solver_type);
    }
  }

  if (section) {
    this->initDOFManager(*section, solver_type);
  } else {
    this->initDOFManager(solver_type);
  }
}

/* -------------------------------------------------------------------------- */
void ModelSolver::initDOFManager(const ID & solver_type) {
  if (solver_type == "explicit") {
    ID id = this->parent_id + ":dof_manager_default";
    this->dof_manager = new DOFManagerDefault(mesh, id, this->parent_memory_id);
  } else if (solver_type == "petsc") {
#if defined(AKANTU_USE_PETSC)
    ID id = this->parent_id + ":dof_manager_petsc";
    this->dof_manager = new DOFManagerPETSc(mesh, id, this->parent_memory_id);
#else
    AKANTU_EXCEPTION(
        "To use PETSc you have to activate it in the compilations options");
#endif
  } else if (solver_type == "mumps") {
#if defined(AKANTU_USE_MUMPS)
    ID id = this->parent_id + ":dof_manager_default";
    this->dof_manager = new DOFManagerDefault(mesh, id, this->parent_memory_id);
#else
    AKANTU_EXCEPTION(
        "To use MUMPS you have to activate it in the compilations options");
#endif
  } else {
    AKANTU_EXCEPTION(
        "To use the solver "
        << solver_type
        << " you will have to code it. This is an unknown solver type.");
  }

    this->setDOFManager(*this->dof_manager);
}

/* -------------------------------------------------------------------------- */
template <typename T> static T getOptionToType(const std::string & opt_str) {
  std::stringstream sstr(opt_str);
  T opt;
  sstr >> opt;

  return opt;
}

/* -------------------------------------------------------------------------- */
void ModelSolver::initDOFManager(const ParserSection & section,
                                 const ID & solver_type) {
  this->initDOFManager(solver_type);
  std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
    sub_sect = section.getSubSections(_st_time_step_solver);

  Parser::const_section_iterator it;
  for (it = sub_sect.first; it != sub_sect.second; ++it) {
    ID solver_id = it->getName();

    std::string str = it->getOption();
    TimeStepSolverType tss_type =
      it->getParameter("type", this->getDefaultSolverType());
    ModelSolverOptions tss_options = this->getDefaultSolverOptions(tss_type);

    std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
      sub_solvers_sect = it->getSubSections(_st_non_linear_solver);
    Parser::const_section_iterator sub_it;
    UInt nb_non_linear_solver_section =
      it->getNbSubSections(_st_non_linear_solver);

    NonLinearSolverType nls_type = tss_options.non_linear_solver_type;

    if (nb_non_linear_solver_section == 1) {
      const ParserSection & nls_section = *(sub_solvers_sect.first);
      nls_type = getOptionToType<NonLinearSolverType>(nls_section.getName());
    } else if (nb_non_linear_solver_section > 0) {
      AKANTU_EXCEPTION("More than one non linear solver are provided for the "
                       "time step solver "
                       << solver_id);
    }

    this->getNewSolver(solver_id, tss_type, nls_type);
    if (nb_non_linear_solver_section == 1) {
      const ParserSection & nls_section = *(sub_solvers_sect.first);
      this->dof_manager->getNonLinearSolver(solver_id).parseSection(
          nls_section);
    }

    std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
      sub_integrator_sect = it->getSubSections(_st_integration_scheme);

    for (sub_it = sub_integrator_sect.first;
         sub_it != sub_integrator_sect.second; ++sub_it) {
      const ParserSection & is_section = *sub_it;
      const ID & dof_id = is_section.getName();

      IntegrationSchemeType it_type = is_section.getParameter(
          "type", tss_options.integration_scheme_type[dof_id]);

      IntegrationScheme::SolutionType s_type = is_section.getParameter(
          "solution_type", tss_options.solution_type[dof_id]);
      this->setIntegrationScheme(solver_id, dof_id, it_type, s_type);
    }

    std::map<ID, IntegrationSchemeType>::const_iterator it =
      tss_options.integration_scheme_type.begin();
    std::map<ID, IntegrationSchemeType>::const_iterator end =
      tss_options.integration_scheme_type.end();
    for (; it != end; ++it) {
      if (!this->hasIntegrationScheme(solver_id, it->first)) {
        this->setIntegrationScheme(solver_id, it->first, it->second,
                                   tss_options.solution_type[it->first]);
      }
    }
  }

  if (section.hasParameter("default_solver")) {
    ID default_solver = section.getParameter("default_solver");
    this->setDefaultSolver(default_solver);
  }
}

/* -------------------------------------------------------------------------- */
TimeStepSolver & ModelSolver::getSolver(const ID & solver_id) {
  ID tmp_solver_id = solver_id;
  if (tmp_solver_id == "")
    tmp_solver_id = this->default_solver_id;

  TimeStepSolver & tss = this->dof_manager->getTimeStepSolver(tmp_solver_id);
  return tss;
}

/* -------------------------------------------------------------------------- */
const TimeStepSolver & ModelSolver::getSolver(const ID & solver_id) const {
  ID tmp_solver_id = solver_id;
  if (solver_id == "")
    tmp_solver_id = this->default_solver_id;

  const TimeStepSolver & tss =
      this->dof_manager->getTimeStepSolver(tmp_solver_id);
  return tss;
}

/* -------------------------------------------------------------------------- */
TimeStepSolver & ModelSolver::getTimeStepSolver(const ID & solver_id) {
  return this->getSolver(solver_id);
}

/* -------------------------------------------------------------------------- */
const TimeStepSolver &
ModelSolver::getTimeStepSolver(const ID & solver_id) const {
  return this->getSolver(solver_id);
}

/* -------------------------------------------------------------------------- */
NonLinearSolver & ModelSolver::getNonLinearSolver(const ID & solver_id) {
  return this->getSolver(solver_id).getNonLinearSolver();
}
/* -------------------------------------------------------------------------- */
const NonLinearSolver &
ModelSolver::getNonLinearSolver(const ID & solver_id) const {
  return this->getSolver(solver_id).getNonLinearSolver();
}

/* -------------------------------------------------------------------------- */
bool ModelSolver::hasSolver(const ID & solver_id) const {
  ID tmp_solver_id = solver_id;
  if (solver_id == "")
    tmp_solver_id = this->default_solver_id;

  return this->dof_manager->hasTimeStepSolver(tmp_solver_id);
}

/* -------------------------------------------------------------------------- */
void ModelSolver::setDefaultSolver(const ID & solver_id) {
  AKANTU_DEBUG_ASSERT(
      this->hasSolver(solver_id),
      "Cannot set the default solver to a solver that does not exists");
  this->default_solver_id = solver_id;
}

/* -------------------------------------------------------------------------- */
void ModelSolver::solveStep(const ID & solver_id) {
  AKANTU_DEBUG_IN();

  TimeStepSolver & tss = this->getSolver(solver_id);
  // make one non linear solve
  tss.solveStep(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ModelSolver::getNewSolver(const ID & solver_id,
                               TimeStepSolverType time_step_solver_type,
                               NonLinearSolverType non_linear_solver_type) {
  if (this->default_solver_id == "") {
    this->default_solver_id = solver_id;
  }

  if (non_linear_solver_type == _nls_auto) {
    switch (time_step_solver_type) {
    case _tsst_dynamic:
    case _tsst_static:
      non_linear_solver_type = _nls_newton_raphson;
      break;
    case _tsst_dynamic_lumped:
      non_linear_solver_type = _nls_lumped;
      break;
    }
  }

  this->initSolver(time_step_solver_type, non_linear_solver_type);

  NonLinearSolver & nls = this->dof_manager->getNewNonLinearSolver(
      solver_id, non_linear_solver_type);

  this->dof_manager->getNewTimeStepSolver(solver_id, time_step_solver_type,
                                          nls);
}

/* -------------------------------------------------------------------------- */
Real ModelSolver::getTimeStep(const ID & solver_id) const {
  const TimeStepSolver & tss = this->getSolver(solver_id);

  return tss.getTimeStep();
}

/* -------------------------------------------------------------------------- */
void ModelSolver::setTimeStep(Real time_step, const ID & solver_id) {
  TimeStepSolver & tss = this->getSolver(solver_id);

  return tss.setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
void ModelSolver::setIntegrationScheme(
    const ID & solver_id, const ID & dof_id,
    const IntegrationSchemeType & integration_scheme_type,
    IntegrationScheme::SolutionType solution_type) {
  TimeStepSolver & tss = this->dof_manager->getTimeStepSolver(solver_id);

  tss.setIntegrationScheme(dof_id, integration_scheme_type, solution_type);
}

/* -------------------------------------------------------------------------- */
bool ModelSolver::hasDefaultSolver() const {
  return (this->default_solver_id != "");
}

/* -------------------------------------------------------------------------- */
bool ModelSolver::hasIntegrationScheme(const ID & solver_id,
                                       const ID & dof_id) const {
  TimeStepSolver & tss = this->dof_manager->getTimeStepSolver(solver_id);
  return tss.hasIntegrationScheme(dof_id);
}

/* -------------------------------------------------------------------------- */
void ModelSolver::predictor() {}

/* -------------------------------------------------------------------------- */
void ModelSolver::corrector() {}

/* -------------------------------------------------------------------------- */
TimeStepSolverType ModelSolver::getDefaultSolverType() const {
  return _tsst_dynamic_lumped;
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions
ModelSolver::getDefaultSolverOptions(__attribute__((unused))
                                     const TimeStepSolverType & type) const {
  ModelSolverOptions options;
  options.non_linear_solver_type = _nls_auto;
  return options;
}

/* -------------------------------------------------------------------------- */
// void ModelSolver::setIntegrationScheme(
//     const ID & solver_id, const ID & dof_id,
//     const IntegrationSchemeType & integration_scheme_type) {}

/* -------------------------------------------------------------------------- */

// /* --------------------------------------------------------------------------
// */
// template <NewmarkBeta::IntegrationSchemeCorrectorType type>
// void SolidMechanicsModel::solve(Array<Real> &increment, Real block_val,
//                                 bool need_factorize, bool
//                                 has_profile_changed) {

//   if(has_profile_changed) {
//     this->initJacobianMatrix();
//     need_factorize = true;
//   }

//   updateResidualInternal(); //doesn't do anything for static

//   if(need_factorize) {
//     Real c = 0.,d = 0.,e = 0.;

//     if(method == _static) {
//       AKANTU_DEBUG_INFO("Solving K inc = r");
//       e = 1.;
//     } else {
//       AKANTU_DEBUG_INFO("Solving (c M + d C + e K) inc = r");

//       NewmarkBeta * nmb_int = dynamic_cast<NewmarkBeta *>(integrator);
//       c = nmb_int->getAccelerationCoefficient<type>(time_step);
//       d = nmb_int->getVelocityCoefficient<type>(time_step);
//       e = nmb_int->getDisplacementCoefficient<type>(time_step);
//     }

//     jacobian_matrix->clear();
//     // J = c M + d C + e K
//     if(stiffness_matrix)
//       jacobian_matrix->add(*stiffness_matrix, e);

//     if(mass_matrix)
//       jacobian_matrix->add(*mass_matrix, c);

// #if !defined(AKANTU_NDEBUG)
//     if(mass_matrix && AKANTU_DEBUG_TEST(dblDump))
//       mass_matrix->saveMatrix("M.mtx");
// #endif

//     if(velocity_damping_matrix)
//       jacobian_matrix->add(*velocity_damping_matrix, d);

//     jacobian_matrix->applyBoundary(*blocked_dofs, block_val);

// #if !defined(AKANTU_NDEBUG)
//     if(AKANTU_DEBUG_TEST(dblDump))
//       jacobian_matrix->saveMatrix("J.mtx");
// #endif
//     solver->factorize();
//   }

//   // if (rhs.getSize() != 0)
//   //  solver->setRHS(rhs);
//   // else

//   solver->setOperators();

//   solver->setRHS(*residual);

//   // solve @f[ J \delta w = r @f]
//   solver->solve(increment);

//   UInt nb_nodes = displacement-> getSize();
//   UInt nb_degree_of_freedom = displacement->getNbComponent() * nb_nodes;

//   bool * blocked_dofs_val = blocked_dofs->storage();
//   Real * increment_val = increment.storage();

//   for (UInt j = 0; j < nb_degree_of_freedom;
//        ++j,++increment_val, ++blocked_dofs_val) {
//     if ((*blocked_dofs_val))
//       *increment_val = 0.0;
//     }

// }

// /* --------------------------------------------------------------------------
// */
// template<SolveConvergenceMethod cmethod, SolveConvergenceCriteria criteria>
// bool SolidMechanicsModel::solveStatic(Real tolerance, UInt max_iteration,
//                                       bool do_not_factorize) {

//   AKANTU_DEBUG_INFO("Solving Ku = f");
//   AKANTU_DEBUG_ASSERT(stiffness_matrix != NULL,
//                       "You should first initialize the implicit solver and
//                       assemble the stiffness matrix by calling
//                       initImplicit");

//   AnalysisMethod analysis_method=method;
//   Real error = 0.;
//   method=_static;
//   bool converged = this->template solveStep<cmethod, criteria>(tolerance,
//   error, max_iteration, do_not_factorize);
//   method=analysis_method;
//   return converged;

// }

// /* --------------------------------------------------------------------------
// */
// template<SolveConvergenceMethod cmethod, SolveConvergenceCriteria criteria>
// bool SolidMechanicsModel::solveStep(Real tolerance,
//                                     UInt max_iteration) {
//   Real error = 0.;
//   return this->template solveStep<cmethod,criteria>(tolerance,
//                                                     error,
//                                                     max_iteration);
// }

// /* --------------------------------------------------------------------------
// */
// template<SolveConvergenceMethod cmethod, SolveConvergenceCriteria criteria>
// bool SolidMechanicsModel::solveStep(Real tolerance, Real & error, UInt
// max_iteration,
//                                     bool do_not_factorize) {
//   EventManager::sendEvent(SolidMechanicsModelEvent::BeforeSolveStepEvent(method));
//   this->implicitPred();
//   this->updateResidual();

//   AKANTU_DEBUG_ASSERT(stiffness_matrix != NULL,
//                       "You should first initialize the implicit solver and
//                       assemble the stiffness matrix");

//   bool need_factorize = !do_not_factorize;

//   if (method==_implicit_dynamic) {
//     AKANTU_DEBUG_ASSERT(mass_matrix != NULL,
//                         "You should first initialize the implicit solver and
//                         assemble the mass matrix");
//   }

//   switch (cmethod) {
//   case _scm_newton_raphson_tangent:
//   case _scm_newton_raphson_tangent_not_computed:
//     break;
//   case _scm_newton_raphson_tangent_modified:
//     this->assembleStiffnessMatrix();
//     break;
//   default:
//     AKANTU_DEBUG_ERROR("The resolution method " << cmethod << " has not been
//     implemented!");
//   }

//   this->n_iter = 0;
//   bool converged = false;
//   error = 0.;
//   if(criteria == _scc_residual) {
//     converged = this->testConvergence<criteria> (tolerance, error);
//     if(converged) return converged;
//   }

//   do {
//     if (cmethod == _scm_newton_raphson_tangent)
//       this->assembleStiffnessMatrix();

//     solve<NewmarkBeta::_displacement_corrector> (*increment, 1.,
//     need_factorize);

//     this->implicitCorr();

//     if(criteria == _scc_residual) this->updateResidual();

//     converged = this->testConvergence<criteria> (tolerance, error);

//     if(criteria == _scc_increment && !converged) this->updateResidual();
//     //this->dump();

//     this->n_iter++;
//     AKANTU_DEBUG_INFO("[" << criteria << "] Convergence iteration "
//                       << std::setw(std::log10(max_iteration)) << this->n_iter
//                       << ": error " << error << (converged ? " < " : " > ")
//                       << tolerance);

//     switch (cmethod) {
//     case _scm_newton_raphson_tangent:
//       need_factorize = true;
//       break;
//     case _scm_newton_raphson_tangent_not_computed:
//     case _scm_newton_raphson_tangent_modified:
//       need_factorize = false;
//       break;
//     default:
//       AKANTU_DEBUG_ERROR("The resolution method " << cmethod << " has not
//       been implemented!");
//     }

//   } while (!converged && this->n_iter < max_iteration);

//   // this makes sure that you have correct strains and stresses after the
//   solveStep function (e.g., for dumping)
//   if(criteria == _scc_increment) this->updateResidual();

//   if (converged) {
//     EventManager::sendEvent(SolidMechanicsModelEvent::AfterSolveStepEvent(method));
//   } else if(this->n_iter == max_iteration) {
//     AKANTU_DEBUG_WARNING("[" << criteria << "] Convergence not reached after
//     "
//                          << std::setw(std::log10(max_iteration)) <<
//                          this->n_iter <<
//                          " iteration" << (this->n_iter == 1 ? "" : "s") <<
//                          "!" << std::endl);
//   }

//   return converged;
// }

// void SolidMechanicsModel::updateResidualInternal() {
//   AKANTU_DEBUG_IN();

//   AKANTU_DEBUG_INFO("Update the residual");
//   // f = f_ext - f_int - Ma - Cv = r - Ma - Cv;

//   if(method != _static) {
//     // f -= Ma
//     if(mass_matrix) {
//       // if full mass_matrix
//       Array<Real> * Ma = new Array<Real>(*acceleration, true, "Ma");
//       *Ma *= *mass_matrix;
//       /// \todo check unit conversion for implicit dynamics
//       //      *Ma /= f_m2a
//       *residual -= *Ma;
//       delete Ma;
//     } else if (mass) {

//       // else lumped mass
//       UInt nb_nodes = acceleration->getSize();
//       UInt nb_degree_of_freedom = acceleration->getNbComponent();

//       Real * mass_val     = mass->storage();
//       Real * accel_val    = acceleration->storage();
//       Real * res_val      = residual->storage();
//       bool * blocked_dofs_val = blocked_dofs->storage();

//       for (UInt n = 0; n < nb_nodes * nb_degree_of_freedom; ++n) {
// 	if(!(*blocked_dofs_val)) {
// 	  *res_val -= *accel_val * *mass_val /f_m2a;
// 	}
// 	blocked_dofs_val++;
// 	res_val++;
// 	mass_val++;
// 	accel_val++;
//       }
//     } else {
//       AKANTU_DEBUG_ERROR("No function called to assemble the mass matrix.");
//     }

//     // f -= Cv
//     if(velocity_damping_matrix) {
//       Array<Real> * Cv = new Array<Real>(*velocity);
//       *Cv *= *velocity_damping_matrix;
//       *residual -= *Cv;
//       delete Cv;
//     }
//   }

//   AKANTU_DEBUG_OUT();
// }

// /* --------------------------------------------------------------------------
// */
// void SolidMechanicsModel::solveLumped(Array<Real> & x,
// 				      const Array<Real> & A,
// 				      const Array<Real> & b,
// 				      const Array<bool> & blocked_dofs,
// 				      Real alpha) {

//   Real * A_val = A.storage();
//   Real * b_val = b.storage();
//   Real * x_val = x.storage();
//   bool * blocked_dofs_val = blocked_dofs.storage();

//   UInt nb_degrees_of_freedom = x.getSize() * x.getNbComponent();

//   for (UInt n = 0; n < nb_degrees_of_freedom; ++n) {
//     if(!(*blocked_dofs_val)) {
//       *x_val = alpha * (*b_val / *A_val);
//     }
//     x_val++;
//     A_val++;
//     b_val++;
//     blocked_dofs_val++;
//   }
// }

/* -------------------------------------------------------------------------- */
// void TimeStepSolverDefault::updateAcceleration() {
//   AKANTU_DEBUG_IN();

//   updateResidualInternal();

//   if (method == _explicit_lumped_mass) {
//     /* residual = residual_{n+1} - M * acceleration_n therefore
//        solution = increment acceleration not acceleration */
//     solveLumped(*increment_acceleration, *mass, *residual, *blocked_dofs,
//                 f_m2a);
//   } else if (method == _explicit_consistent_mass) {
//     solve<NewmarkBeta::_acceleration_corrector>(*increment_acceleration);
//   }

//   AKANTU_DEBUG_OUT();
// }

} // akantu
