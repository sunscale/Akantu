/**
 * @file   model_solver.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Aug 18 2015
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Implementation of ModelSolver
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
template <typename T> static T getOptionToType(const std::string & opt_str) {
  std::stringstream sstr(opt_str);
  T opt;
  sstr >> opt;

  return opt;
}

/* -------------------------------------------------------------------------- */
ModelSolver::ModelSolver(Mesh & mesh, const ModelType & type, const ID & id)
    : Parsable(ParserType::_model, id), model_type(type), parent_id(id),
      mesh(mesh), dof_manager(nullptr) {}

/* -------------------------------------------------------------------------- */
ModelSolver::~ModelSolver() = default;

/* -------------------------------------------------------------------------- */
std::tuple<ParserSection, bool> ModelSolver::getParserSection() {
  auto sub_sections = getStaticParser().getSubSections(ParserType::_model);

  auto it = std::find_if(
      sub_sections.begin(), sub_sections.end(), [&](auto && section) {
        auto type = getOptionToType<ModelType>(section.getName());
        // default id should be the model type if not defined
        std::string name = section.getParameter("name", this->parent_id);
        return type == model_type and name == this->parent_id;
      });

  if (it == sub_sections.end()) {
    return std::make_tuple(ParserSection(), true);
  }

  return std::make_tuple(*it, false);
}

/* -------------------------------------------------------------------------- */
void ModelSolver::initDOFManager() {
  // default without external solver activated at compilation same as mumps that
  // is the historical solver but with only the lumped solver
  ID solver_type = "default";

#if defined(AKANTU_USE_MUMPS)
  solver_type = "default";
#elif defined(AKANTU_USE_PETSC)
  solver_type = "petsc";
#endif

  ParserSection section;
  bool is_empty;
  std::tie(section, is_empty) = this->getParserSection();

  if (not is_empty) {
    solver_type = section.getOption(solver_type);
    this->initDOFManager(section, solver_type);
  } else {
    this->initDOFManager(solver_type);
  }
}

/* -------------------------------------------------------------------------- */
void ModelSolver::initDOFManager(const ID & solver_type) {
  try {
    this->dof_manager = DOFManagerFactory::getInstance().allocate(
        solver_type, mesh, this->parent_id + ":dof_manager_" + solver_type);
  } catch (...) {
    AKANTU_EXCEPTION(
        "To use the solver "
        << solver_type
        << " you will have to code it. This is an unknown solver type.");
  }

  this->setDOFManager(*this->dof_manager);
}

/* -------------------------------------------------------------------------- */
void ModelSolver::initDOFManager(const ParserSection & section,
                                 const ID & solver_type) {
  this->initDOFManager(solver_type);
  auto sub_sections = section.getSubSections(ParserType::_time_step_solver);

  // parsing the time step solvers
  for (auto && section : sub_sections) {
    ID type = section.getName();
    ID solver_id = section.getParameter("name", type);

    auto tss_type = getOptionToType<TimeStepSolverType>(type);
    auto tss_options = this->getDefaultSolverOptions(tss_type);

    auto sub_solvers_sect =
        section.getSubSections(ParserType::_non_linear_solver);
    auto nb_non_linear_solver_section =
        section.getNbSubSections(ParserType::_non_linear_solver);

    auto nls_type = tss_options.non_linear_solver_type;

    if (nb_non_linear_solver_section == 1) {
      auto && nls_section = *(sub_solvers_sect.first);
      nls_type = getOptionToType<NonLinearSolverType>(nls_section.getName());
    } else if (nb_non_linear_solver_section > 0) {
      AKANTU_EXCEPTION("More than one non linear solver are provided for the "
                       "time step solver "
                       << solver_id);
    }

    this->getNewSolver(solver_id, tss_type, nls_type);
    if (nb_non_linear_solver_section == 1) {
      const auto & nls_section = *(sub_solvers_sect.first);
      this->dof_manager->getNonLinearSolver(solver_id).parseSection(
          nls_section);
    }

    auto sub_integrator_sections =
        section.getSubSections(ParserType::_integration_scheme);

    for (auto && is_section : sub_integrator_sections) {
      const auto & dof_type_str = is_section.getName();
      ID dof_id;
      try {
        ID tmp = is_section.getParameter("name");
        dof_id = tmp;
      } catch (...) {
        AKANTU_EXCEPTION("No degree of freedom name specified for the "
                         "integration scheme of type "
                         << dof_type_str);
      }

      auto it_type = getOptionToType<IntegrationSchemeType>(dof_type_str);

      IntegrationScheme::SolutionType s_type = is_section.getParameter(
          "solution_type", tss_options.solution_type[dof_id]);
      this->setIntegrationScheme(solver_id, dof_id, it_type, s_type);
    }

    for (auto & is_type : tss_options.integration_scheme_type) {
      if (!this->hasIntegrationScheme(solver_id, is_type.first)) {
        this->setIntegrationScheme(solver_id, is_type.first, is_type.second,
                                   tss_options.solution_type[is_type.first]);
      }
    }
  }

  if (section.hasParameter("default_solver")) {
    ID default_solver = section.getParameter("default_solver");
    if (this->hasSolver(default_solver)) {
      this->setDefaultSolver(default_solver);
    } else {
      AKANTU_EXCEPTION(
          "The solver \""
          << default_solver
          << "\" was not created, it cannot be set as default solver");
    }
  }
}

/* -------------------------------------------------------------------------- */
TimeStepSolver & ModelSolver::getSolver(const ID & solver_id) {
  ID tmp_solver_id = solver_id;
  if (tmp_solver_id.empty()) {
    tmp_solver_id = this->default_solver_id;
  }

  TimeStepSolver & tss = this->dof_manager->getTimeStepSolver(tmp_solver_id);
  return tss;
}

/* -------------------------------------------------------------------------- */
const TimeStepSolver & ModelSolver::getSolver(const ID & solver_id) const {
  ID tmp_solver_id = solver_id;
  if (solver_id.empty()) {
    tmp_solver_id = this->default_solver_id;
  }

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
  if (solver_id.empty()) {
    tmp_solver_id = this->default_solver_id;
  }

  if (not this->dof_manager) {
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
void ModelSolver::solveStep(SolverCallback & callback, const ID & solver_id) {
  AKANTU_DEBUG_IN();

  TimeStepSolver & tss = this->getSolver(solver_id);
  // make one non linear solve
  tss.solveStep(callback);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ModelSolver::solveStep(const ID & solver_id) {
  solveStep(*this, solver_id);
}

/* -------------------------------------------------------------------------- */
void ModelSolver::getNewSolver(const ID & solver_id,
                               TimeStepSolverType time_step_solver_type,
                               NonLinearSolverType non_linear_solver_type) {
  if (this->default_solver_id.empty()) {
    this->default_solver_id = solver_id;
  }

  if (non_linear_solver_type == NonLinearSolverType::_auto) {
    switch (time_step_solver_type) {
    case TimeStepSolverType::_dynamic:
    case TimeStepSolverType::_static:
      non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      break;
    case TimeStepSolverType::_dynamic_lumped:
      non_linear_solver_type = NonLinearSolverType::_lumped;
      break;
    case TimeStepSolverType::_not_defined:
      AKANTU_EXCEPTION(time_step_solver_type
                       << " is not a valid time step solver type");
      break;
    }
  }

  this->initSolver(time_step_solver_type, non_linear_solver_type);

  NonLinearSolver & nls = this->dof_manager->getNewNonLinearSolver(
      solver_id, non_linear_solver_type);

  this->dof_manager->getNewTimeStepSolver(solver_id, time_step_solver_type, nls,
                                          *this);
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
  return (not this->default_solver_id.empty());
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
  return TimeStepSolverType::_dynamic_lumped;
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions
ModelSolver::getDefaultSolverOptions(__attribute__((unused))
                                     const TimeStepSolverType & type) const {
  ModelSolverOptions options;
  options.non_linear_solver_type = NonLinearSolverType::_auto;
  return options;
}

} // namespace akantu
