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
#include "dof_manager_default.hh"
#include "mesh.hh"
#include "non_linear_solver.hh"
#include "time_step_solver.hh"

#if defined(AKANTU_USE_PETSC)
#include "dof_manager_petsc.hh"
#endif

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ModelSolver::ModelSolver(Mesh & mesh, const ID & id, UInt memory_id)
    : Parsable(_st_model_solver, id), SolverCallback(), parent_id(id),
      parent_memory_id(memory_id), mesh(mesh), dof_manager(nullptr),
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

  const ParserSection * section = nullptr;
  Parser::const_section_iterator it;
  for (it = sub_sect.first; it != sub_sect.second && section == nullptr; ++it) {
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

    // std::string str = it->getOption();
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

  if(not this->dof_manager) {
    AKANTU_EXCEPTION("No DOF manager was initialized");
  }
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
    case _tsst_not_defined:
      AKANTU_EXCEPTION(time_step_solver_type
                       << " is not a valid time step solver type");
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

} // namespace akantu
