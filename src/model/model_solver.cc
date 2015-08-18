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
#if defined(AKANTU_USE_MUMPS)
#include "dof_manager_default.hh"
#endif
#if defined(AKANTU_USE_PETSC)
#include "dof_manager_petsc.hh"
#endif

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

ModelSolver::ModelSolver(const ID & id)
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
    this->dof_manager = new DOFManagerPETSc();
#else
    AKANTU_EXCEPTION(
        "To use PETSc you have to activate it in the compilations options");
#endif
  } else if (solver_type == "mumps") {
#if defined(AKANTU_USE_MUMPS)
    this->dof_manager = new DOFManagerDefault();
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

__END_AKANTU__
