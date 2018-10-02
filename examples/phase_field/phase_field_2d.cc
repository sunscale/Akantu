/**
 * @file   phase_field_static_2d.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Oct 1 2018
 *
 * @brief  test of the class PhaseFieldModel on the 2d square
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "phase_field_model.hh"
#include "solid_mechanics_model.hh"
#include "solid_phase_coupler.hh"
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
const UInt spatial_dimension = 2;
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  initialize("materia.dat", argc, argv);

  // create mesh
  Mesh mesh(spatial_dimension);
  mesh.read("square.msh");
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  // Phase field model initialization
  PhaseFieldModel pfm(mesh);
  pfm.initFull(_analysis_method = _static);

  auto & solver = pfm.getNonLinearSolver();
  solver.set("max_iterations", 1000);
  solver.set("threshold", 1e-8);
  solver.set("convergence_type", _scc_solution);
  
  // solid mechanics model initialization
  SolidMechanicsModel smm(mesh);
  smm.initFull(_analysis_method  = _static);

  smm.applyBC(BC::Dirichlet::FixedValue(0., _y), "bottom");
  smm.applyBC(BC::Dirichlet::FixedValue(0., _x), "bottom");
  smm.applyBC(BC::Dirichlet::FixedValue(0., _x), "left");
  smm.applyBC(BC::Dirichlet::FixedValue(0., _x), "right");

  smm.setBasename(        "square");
  smm.addDumpFieldVector( "displacement");
  smm.addDumpFieldVector( "residual");
  smm.addDumpField(       "stress");
  smm.addDumpField(       "grad_u");
  smm.addDumpField(       "damage");
  smm.addDumpField(       "blocked_dofs");
  smm.dump();

  solver = smm.getNonLinearSolver();
  solver.set("max_iterations", 1000);
  solver.set("threshold", 1e-8);
  solver.set("convergence_type", _scc_solution);
  
  // coupling of models
  SolidPhaseCoupler coupler<smm, pfm>();
  // assign the type of solver for coupler whether it is staggered or
  // monolithic
  // auto & scheme = coupler.getCouplingScheme();
  // scheme.set("max_iterations", 1000);
  // scheme.set("threshold", 1e-8);
  // scheme.set("convergence_dof", _damage);

  UInt nbSteps   = 1000;
  Real increment = 1.e-4;
  
  for (UInt s = 1; s < nbSteps; ++s) {
    smm.applyBC(BC::Dirichlet::IncrementValue(increment, _y), "top");
    coupler.solveStep();
    smm.dump();
  }

  finialize;
  return EXIT_SUCCESS;
}

