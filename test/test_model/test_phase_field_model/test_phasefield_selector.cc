/**
 * @file   test_phasefield_selector.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Jun 18 2020
 * @date last modification: Thu Jun 18 2020
 *
 * @brief  Test for phasefield selector
 *
 * @section LICENSE
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

#include "aka_common.hh"
#include "phase_field_model.hh"
#include "phasefield_exponential.hh"

using namespace akantu;

int main(int argc, char *argv[]){

  initialize("phasefield_selector.dat", argc, argv);

  Math::setTolerance(1e-8);

  Mesh mesh(1);
  mesh.read("phasefield_selector.msh");

  PhaseFieldModel model(mesh);
  auto && selector = std::make_shared<MeshDataPhaseFieldSelector<std::string>>(
									       "physical_names", model);
  model.setPhaseFieldSelector(selector);

  model.initFull();

  PhaseField & chocolate = model.getPhaseField("chocolate");
  PhaseField & chewing_gum = model.getPhaseField("chewing-gum");
  PhaseField & candy = model.getPhaseField("candy");

  const Real gc = chocolate.getParam("gc");
  

  UInt chocolate_element = chocolate.getElementFilter(_segment_2)(0, 0);
  UInt chewing_gum_element = chewing_gum.getElementFilter(_segment_2)(0, 0);
  UInt candy_element = candy.getElementFilter(_segment_2)(0, 0);
  
  if (chocolate_element != 0 || chewing_gum_element != 1 || candy_element != 2 ) {
    return EXIT_FAILURE;
  }
 
  return EXIT_SUCCESS;
}


