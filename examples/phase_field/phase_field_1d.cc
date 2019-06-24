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
#include "non_linear_solver.hh"
#include "phase_field_model.hh"
#include "solid_mechanics_model.hh"
#include "solid_phase_coupler.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
const UInt spatial_dimension = 1;
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  initialize("material_1d.dat", argc, argv);

  const Real max_displacement = 0.05;

  const UInt nbSteps = 50 ;
  const UInt unload_start = 200; // add to be evenmultiple of 10
  Real time_step = max_displacement/nbSteps;
  
  Mesh mesh(spatial_dimension);
  mesh.read("1d_1elem_bar.msh");
 
  PhaseFieldModel phase(mesh);
  phase.initFull(_analysis_method = _static);

  // solid mechanics model initialization
  SolidMechanicsModel solid(mesh);
  solid.initFull(_analysis_method  = _static);

  solid.applyBC(BC::Dirichlet::FixedValue(0., _x), "Left");
  solid.applyBC(BC::Dirichlet::FixedValue(0., _x), "Right"); 

  solid.setBaseName(        "square");
  solid.addDumpFieldVector( "displacement");
  solid.addDumpFieldVector( "internal_force");
  solid.addDumpField(       "stress");
  solid.addDumpField(       "grad_u");
  solid.addDumpField(       "damage");
  solid.addDumpField(       "blocked_dofs");
  solid.dump();

  auto & solid_solver = solid.getNonLinearSolver();
  solid_solver.set("max_iterations", 1000);
  solid_solver.set("threshold", 1e-8);
  solid_solver.set("convergence_type", SolveConvergenceCriteria::_residual);
  
  // coupling of models
  SolidPhaseCoupler<SolidMechanicsModel, PhaseFieldModel> coupler(solid, phase);
  
  Real stress_homogeneous;
  Real young_unload;
  Real gc = 0.00014e-2;
  Real l0 = 1./8;
  Real Young = 1.0;
  
  for (UInt s = 1; s < nbSteps; ++s) {
    //Increasing loading
    if( s<unload_start || s>(2.1*unload_start)){
      solid.applyBC(BC::Dirichlet::IncrementValue(-time_step,_x),"Left");
      solid.applyBC(BC::Dirichlet::IncrementValue(time_step,_x),"Right");
    }
    else{
      solid.applyBC(BC::Dirichlet::IncrementValue(time_step,_x),"Left");
      solid.applyBC(BC::Dirichlet::IncrementValue(-time_step,_x),"Right");
    }

    coupler.solve();
    
    Array<Real> & stress = solid.getMaterial("solid").getArray<Real>("stress", _segment_2);
    Array<Real> & damage = solid.getMaterial("solid").getArray<Real>("damage", _segment_2);
    Array<Real> & grad_u = solid.getMaterial("solid").getArray<Real>("grad_u", _segment_2);
	
    solid.dump();

    if(s==(unload_start-1))
      young_unload = stress(0,0)/grad_u(0,0);
    if( s<unload_start || s>(3.2*unload_start)) {
      // verification that the stress match the 1D analytical homogeneous stress
      // from Borden et al. CMAME vol. 217-220, page 77-95 (2012)
      stress_homogeneous = pow(l0/gc*Young*pow(grad_u(0,0),2)+1,-2)*Young*grad_u(0,0);
      std::cout << stress_homogeneous << std::endl;
      std::cout << stress(0,0) << " ---- "  << damage(0,0) << std::endl;
          
      //if( (std::abs(stress_homogeneous-stress(0,0))/stress_homogeneous) > 1e-9)
      //	return EXIT_FAILURE;
    }
    else {
      Real sig_outof_eps = std::abs(stress(0,0)/grad_u(0,0));
      if( (std::abs(grad_u(0,0)) > 1e-9)
	  && ( (grad_u(0,0)>0 && (std::abs(sig_outof_eps-young_unload)/young_unload) > 1e-9)
	       || ( grad_u(0,0)<0 && std::abs(sig_outof_eps-Young)/Young > 1e-9 ) ) ){
	std::cout << s << "," <<grad_u(0,0) << "," << stress(0,0) / grad_u(0,0) << "," << Young << "," << std::abs(sig_outof_eps-young_unload)/young_unload << "," << std::abs(sig_outof_eps-Young)/Young << std::endl;
	return EXIT_FAILURE;
       }
    }				      
								       
    std::cout << "Step " << s << "/" << nbSteps << std::endl;
  }

  finalize();
  return EXIT_SUCCESS;
}

