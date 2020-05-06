/**
 * @file   test_material_selector.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri May 01 2015
 * @date last modification: Tue Dec 05 2017
 *
 * @brief  Test for material selector
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

#include "aka_common.hh"
#include "solid_mechanics_model.hh"

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material_selector.dat", argc, argv);

  Math::setTolerance(1e-8);

  Mesh mesh(1);
  mesh.read("material_selector.msh");

  SolidMechanicsModel model(mesh);
  auto && selector = std::make_shared<MeshDataMaterialSelector<std::string>>(
      "physical_names", model);
  model.setMaterialSelector(selector);

  model.initFull();

  Material & chocolate = model.getMaterial("chocolate");
  Material & chewing_gum = model.getMaterial("chewing-gum");
  Material & candy = model.getMaterial("candy");

  UInt chocolate_element = chocolate.getElementFilter(_segment_2)(0, 0);
  UInt chewing_gum_element = chewing_gum.getElementFilter(_segment_2)(0, 0);
  UInt candy_element = candy.getElementFilter(_segment_2)(0, 0);

  if (chocolate_element != 0 || chewing_gum_element != 1 || candy_element != 2)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
