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

#if defined(AKANTU_USE_MPI)
#include "mpi_type_wrapper.hh"
#endif
#if defined(AKANTU_USE_MUMPS)
#include "dof_manager_default.hh"
#endif
#if defined(AKANTU_USE_PETSC)
#include "dof_manager_petsc.hh"
#endif

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

ModelSolver::ModelSolver(const Mesh & mesh, const ID & id, UInt memory_id)
    : Parsable(_st_solver, id), dof_manager(NULL) {
  std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
      sub_sect = getStaticParser().getSubSections(_st_solver);

  if (sub_sect.first != sub_sect.second) {
    AKANTU_EXCEPTION("More than on solver section present in the input file");
  }

  const ParserSection & section = *sub_sect.first;
  std::string solver_type = section.getName();

  if (solver_type == "petsc") {
#if defined(AKANTU_USE_PETSC)
    this->dof_manager = new DOFManagerPETSc(mesh, id + ":dof_manager_petsc", memory_id);
#else
    AKANTU_EXCEPTION(
        "To use PETSc you have to activate it in the compilations options");
#endif
  } else if (solver_type == "mumps") {
#if defined(AKANTU_USE_MUMPS)
    this->dof_manager =
        new DOFManagerDefault(mesh, id + ":dof_manager_default", memory_id);
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
}

/* -------------------------------------------------------------------------- */
void ModelSolver::initTimeStepSolver(
    const ID & time_step_solver_id, const ID & dof_id,
    const TimeStepSolverType & time_step_solver_type) {
  if(this->default_time_step_solver == "") {
    this->default_time_step_solver = time_step_solver_id;
  }

  this->dof_manager->getNewTimeStepSolver(time_step_solver_id, dof_id,
                                          time_step_solver_type);
}

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

__END_AKANTU__
