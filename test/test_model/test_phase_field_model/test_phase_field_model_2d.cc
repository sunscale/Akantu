/**
 * @file   tets_phase_field_2d.cc
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
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
const UInt spatial_dimension = 2;
/* -------------------------------------------------------------------------- */
void applyStrainOnQuadPoints(PhaseFieldModel &);
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]) {

  initialize("material.dat", argc, argv);
  
  Mesh mesh(spatial_dimension);
  mesh.read("test_one_element.msh");

  PhaseFieldModel model(mesh);
  model.initFull(_analysis_method = _static);

  model.setBaseName("one_element");
  model.addDumpFieldVector("damage");
  model.dump();

  UInt nbSteps = 1000;

  for (UInt s = 0; s < nbSteps; ++s) {
    applyStrainOnQuadPoints(model);
    model.solveStep();
    model.dump();
  }
  
  finalize();

  return EXIT_SUCCESS;

}


/* -------------------------------------------------------------------------- */
void applyStrainOnQuadPoints(PhaseFieldModel & model) {
  auto & strain_on_qpoints = model.getStrain();

  auto & mesh = model.getMesh();
  auto & fem = model.getFEEngine(); 
  
  for (auto & type: mesh.elementTypes(spatial_dimension, _not_ghost) ) {
    auto & strain_on_qpoints_vect = strain_on_qpoints(type, _not_ghost);
    UInt nb_quad_points = fem.getNbIntegrationPoints(type, _not_ghost);
    UInt nb_elements = mesh.getNbElement(type);
    
    Array<Real> quad_coords(nb_elements * nb_quad_points, spatial_dimension);
    fem.interpolateOnIntegrationPoints(mesh.getNodes(), quad_coords, spatial_dimension, type, _not_ghost);

    for (auto && values:
	   zip(make_view(quad_coords), 
	       make_view(strain_on_qpoints_vect, spatial_dimension, spatial_dimension))) {
      auto & coord  = std::get<0>(values);
      auto & strain = std::get<1>(values);

      strain(1,1) += 1.e-4;
    }
  }
}
