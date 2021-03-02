/**
 * @file   test_phase_field_explicit.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Feb 28 2021
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
#include "aka_common.hh"
#include "non_linear_solver.hh"
#include "coupler_solid_phasefield.hh"
#include "solid_mechanics_model.hh"
#include "phase_field_model.hh"
#include "material.hh"
#include "material_phasefield.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <fstream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
const UInt spatial_dimension = 2;

/* -------------------------------------------------------------------------- */
void applyDisplacement(SolidMechanicsModel &, Real &);
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]) {

  std::ofstream os("data-explicit.csv");
  os << "#strain stress damage analytical_sigma analytical_damage error_stress error_damage" << std::endl;
  
  initialize("material_coupling.dat", argc, argv);
  
  Mesh mesh(spatial_dimension);
  mesh.read("test_one_element.msh");

  CouplerSolidPhaseField coupler(mesh);
  auto & model = coupler.getSolidMechanicsModel();
  auto & phase = coupler.getPhaseFieldModel();
  
  model.initFull(_analysis_method = _explicit_lumped_mass);

  Real time_factor = 0.8;
  Real stable_time_step = model.getStableTimeStep() * time_factor;
  model.setTimeStep(stable_time_step);

  auto && selector = std::make_shared<MeshDataPhaseFieldSelector<std::string>>(
	  "physical_names", phase);
  phase.setPhaseFieldSelector(selector);
  phase.initFull(_analysis_method = _static);
  
  model.setBaseName("phase_solid");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpFieldVector("displacement");
  model.addDumpField("damage");
  model.dump();
  
  UInt nbSteps = 1000;
  Real increment = 1e-4;

  auto & stress = model.getMaterial(0).getArray<Real>("stress", _quadrangle_4);
  auto & damage = model.getMaterial(0).getArray<Real>("damage", _quadrangle_4);

  Real analytical_damage{0.};
  Real analytical_sigma{0.};

  auto & phasefield = phase.getPhaseField(0);
    
  const Real E   = phasefield.getParam("E");
  const Real nu  = phasefield.getParam("nu");
  Real c22 = E*(1-nu)/((1+nu)*(1-2*nu));

  const Real gc = phasefield.getParam("gc");
  const Real l0 = phasefield.getParam("l0");

  Real error_stress{0.};
  Real error_damage{0.};
 
  for (UInt s = 0; s < nbSteps; ++s) {
    Real axial_strain = increment * s;
    applyDisplacement(model, axial_strain);

    coupler.solve("explicit_lumped", "static");
    
    analytical_damage = axial_strain*axial_strain*c22/(gc/l0 + axial_strain*axial_strain*c22);
    analytical_sigma  = c22*axial_strain*(1-analytical_damage)*(1-analytical_damage);

    error_stress = std::abs(analytical_sigma - stress(0, 3))/analytical_sigma;
    
    error_damage = std::abs(analytical_damage - damage(0))/analytical_damage;

    if (error_damage > 0.01) {
      return EXIT_FAILURE;
    }
    
    os << axial_strain << " " << stress(0, 3) << " " << damage(0) << " "
       << analytical_sigma << " " << analytical_damage << " "  <<
      error_stress  << " " << error_damage << std::endl;

    model.dump();
  }

  os.close();
  finalize();

 

  return EXIT_SUCCESS;

}


/* -------------------------------------------------------------------------- */
void applyDisplacement(SolidMechanicsModel & model, Real & increment) {
  auto & displacement = model.getDisplacement();

  auto & positions = model.getMesh().getNodes();
  auto & blocked_dofs = model.getBlockedDOFs();

  
  for (UInt n = 0; n < model.getMesh().getNbNodes(); ++n) {
    if (positions(n, 1) == -0.5) {
      displacement(n, 0) = 0;
      displacement(n, 1) = 0;
      blocked_dofs(n, 0) = true;
      blocked_dofs(n ,1) = true;
    }
    else {
      displacement(n, 0) = 0;
      displacement(n, 1) = increment;
      blocked_dofs(n, 0) = true;
      blocked_dofs(n ,1) = true;
    }
  }
}
