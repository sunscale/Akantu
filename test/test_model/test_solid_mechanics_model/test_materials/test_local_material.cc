/**
 * @file   test_local_material.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Clement Roux <clement.roux@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Thu Dec 14 2017
 *
 * @brief  test of the class SolidMechanicsModel with custom local damage on a
 * notched plate
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <iostream>
/* -------------------------------------------------------------------------- */
#include "local_material_damage.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  akantu::initialize("material.dat", argc, argv);
  UInt max_steps = 1100;

  const UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("mesh_section_gap.msh");

  /// model initialization
  MaterialFactory::getInstance().registerAllocator(
      "local_damage",
      [](UInt, const ID &, SolidMechanicsModel & model,
         const ID & id) -> std::unique_ptr<Material> {
        return std::make_unique<LocalMaterialDamage>(model, id);
      });

  SolidMechanicsModel model(mesh);
  model.initFull();

  std::cout << model.getMaterial(0) << std::endl;

  model.addDumpField("damage");

  model.addDumpField("strain");
  model.addDumpField("stress");
  model.addDumpFieldVector("displacement");
  model.addDumpFieldVector("external_force");
  model.addDumpFieldVector("internal_force");
  model.dump();

  Real time_step = model.getStableTimeStep();
  model.setTimeStep(time_step / 2.5);

  /// Dirichlet boundary conditions
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "Fixed");
  // model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "Fixed");
  Matrix<Real> stress(2, 2);
  stress.eye(5e7);
  model.applyBC(BC::Neumann::FromHigherDim(stress), "Traction");

  for (UInt s = 0; s < max_steps; ++s)
    model.solveStep();

  model.dump();

  // This should throw a bad_cast if not the proper material
  auto & mat =
      dynamic_cast<LocalMaterialDamage &>(model.getMaterial("concrete"));
  const auto & filter = mat.getElementFilter();
  for (auto & type : filter.elementTypes(spatial_dimension)) {
    std::cout << mat.getDamage(type) << std::endl;
  }

  // This part of the test is to mesh dependent and as nothing to do with the
  // fact that we can create a user defined material or not

  // const auto & lower_bounds = mesh.getLowerBounds();
  // const auto & upper_bounds = mesh.getUpperBounds();
  // Real L = upper_bounds(_x) - lower_bounds(_x);
  // Real H = upper_bounds(_y) - lower_bounds(_y);

  // const auto & filter = model.getMaterial("concrete").getElementFilter();

  // Vector<Real> barycenter(spatial_dimension);

  // for (auto & type : filter.elementTypes(spatial_dimension)) {
  //   UInt nb_elem = mesh.getNbElement(type);
  //   const UInt nb_gp = model.getFEEngine().getNbIntegrationPoints(type);
  //   const auto & material_damage_array =
  //       model.getMaterial(0).getArray<Real>("damage", type);
  //   UInt cpt = 0;
  //   for (auto nel : arange(nb_elem)) {
  //     mesh.getBarycenter({type, nel, _not_ghost}, barycenter);
  //     if ((std::abs(barycenter(_x) - (L / 2) + 0.025) < 0.025) &&
  //         (std::abs(barycenter(_y) - (H / 2) + 0.045) < 0.045)) {
  //       if (material_damage_array(cpt, 0) < 0.9) {
  //         std::terminate();
  //       } else {
  //         std::cout << "element " << nel << " is correctly broken" <<
  //         std::endl;
  //       }
  //     }

  //     cpt += nb_gp;
  //   }
  // }

  akantu::finalize();
  return 0;
}
