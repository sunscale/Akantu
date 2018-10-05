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
  
  PhaseFieldModel pfm(mesh);
  pfm.initFull(_analysis_method = _static);

  if (prank == 0) {
    std::cout << pfm << std::endl;
  }
  
  auto & pfm_solver = pfm.getNonLinearSolver();
  pfm_solver.set("max_iterations", 1000);
  pfm_solver.set("threshold", 1e-3);
  pfm_solver.set("convergence_type", _scc_solution);
  
  SolidMechanicsModel smm(mesh);
  smm.initFull(_analysis_method  = _static);

  smm.applyBC(BC::Dirichlet::FixedValue(0., _y), "bottom");
  smm.applyBC(BC::Dirichlet::FixedValue(0., _x), "bottom");
  smm.applyBC(BC::Dirichlet::FixedValue(0., _x), "left");
  smm.applyBC(BC::Dirichlet::FixedValue(0., _x), "right");

  smm.setBaseName(        "square");
  smm.addDumpFieldVector( "displacement");
  smm.addDumpFieldVector( "internal_force");
  smm.addDumpField(       "stress");
  smm.addDumpField(       "grad_u");
  smm.addDumpField(       "damage");
  smm.addDumpField(       "blocked_dofs");
  smm.dump();

  auto & smm_solver = smm.getNonLinearSolver();
  smm_solver.set("max_iterations", 1000);
  smm_solver.set("threshold", 1e-8);
  smm_solver.set("convergence_type", _scc_solution);
  
  SolidPhaseCoupler<SolidMechanicsModel, PhaseFieldModel> coupler(smm, pfm);

  UInt nbSteps   = 100;
  Real increment = 1.e-4;
  
  for (UInt s = 1; s < nbSteps; ++s) {
    smm.applyBC(BC::Dirichlet::IncrementValue(increment, _y), "top");
    coupler.solve();
    smm.dump();  
    if (prank == 0) {
      std::cout << "Step " << s << "/" << nbSteps << std::endl;
    }
    
  }

  finalize();
  return EXIT_SUCCESS;
}

