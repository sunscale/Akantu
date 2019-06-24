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
#include "communicator.hh"
#include "non_linear_solver.hh"
#include "phase_field_model.hh"
#include "solid_mechanics_model.hh"
#include "solid_phase_coupler.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
const UInt spatial_dimension = 2;
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  initialize("material.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  
  const auto & comm = Communicator::getStaticCommunicator();
  Int prank = comm.whoAmI();

  if (prank == 0) {
    mesh.read("square.msh");
  }

  mesh.distribute();
  
  PhaseFieldModel phase(mesh);
  phase.initFull(_analysis_method = _static);

  if (prank == 0) {
    std::cout << phase << std::endl;
  }
  
  auto & pfm_solver = phase.getNonLinearSolver();
  pfm_solver.set("max_iterations", 1000);
  pfm_solver.set("threshold", 1e-3);
  pfm_solver.set("convergence_type", SolveConvergenceCriteria::_solution);
  
  SolidMechanicsModel solid(mesh);
  solid.initFull(_analysis_method  = _static);

  solid.applyBC(BC::Dirichlet::FixedValue(0., _y), "bottom");
  solid.applyBC(BC::Dirichlet::FixedValue(0., _x), "bottom");
  solid.applyBC(BC::Dirichlet::FixedValue(0., _x), "left");
  solid.applyBC(BC::Dirichlet::FixedValue(0., _x), "right");

  solid.setBaseName(        "square");
  solid.addDumpFieldVector( "displacement");
  solid.addDumpFieldVector( "internal_force");
  solid.addDumpField(       "stress");
  solid.addDumpField(       "grad_u");
  solid.addDumpField(       "damage");
  solid.addDumpField(       "blocked_dofs");
  solid.dump();

  auto & smm_solver = solid.getNonLinearSolver();
  smm_solver.set("max_iterations", 1000);
  smm_solver.set("threshold", 1e-8);
  smm_solver.set("convergence_type", SolveConvergenceCriteria::_solution);
  
  SolidPhaseCoupler<SolidMechanicsModel, PhaseFieldModel> coupler(solid, phase);

  UInt nbSteps   = 100;
  Real increment = 1.e-4;
  
  for (UInt s = 1; s < nbSteps; ++s) {
    solid.applyBC(BC::Dirichlet::IncrementValue(increment, _y), "top");
    coupler.solve();
    solid.dump();  
    if (prank == 0) {
      std::cout << "Step " << s << "/" << nbSteps << std::endl;
    }
    
  }

  finalize();
  return EXIT_SUCCESS;
}

