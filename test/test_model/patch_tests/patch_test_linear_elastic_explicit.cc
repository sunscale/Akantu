/**
 * @file   patch_test_linear_elastic_explicit.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jan 30 2018
 *
 * @brief  patch test solid mechanics explicit
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "patch_test_linear_solid_mechanics_fixture.hh"
/* -------------------------------------------------------------------------- */

TYPED_TEST(TestPatchTestSMMLinear, Explicit) {
  std::string filename = "material_check_stress_plane_stress.dat";
  if (this->plane_strain)
    filename = "material_check_stress_plane_strain.dat";

  this->initModel(_explicit_lumped_mass, filename);

  const auto & coordinates = this->mesh->getNodes();
  auto & displacement = this->model->getDisplacement();
  // set the position of all nodes to the static solution
  for (auto && tuple : zip(make_view(coordinates, this->dim),
                           make_view(displacement, this->dim))) {
    this->setLinearDOF(std::get<1>(tuple), std::get<0>(tuple));
  }
  for (UInt s = 0; s < 100; ++s) {
    this->model->solveStep();
  }

  auto ekin = this->model->getEnergy("kinetic");
  EXPECT_NEAR(0, ekin, 1e-16);

  this->checkAll();

#define debug 0
#if debug
  this->model->setBaseName(std::to_string(this->type) + "_explicit");
  this->model->addDumpField("stress");
  this->model->addDumpField("grad_u");
  this->model->addDumpFieldVector("internal_force");
  this->model->addDumpFieldVector("external_force");
  this->model->addDumpField("blocked_dofs");
  this->model->addDumpFieldVector("displacement");
  this->model->dump();
#endif

}

/* -------------------------------------------------------------------------- */
TYPED_TEST(TestPatchTestSMMLinear, ExplicitFiniteDeformation) {
  std::string filename = "material_check_stress_plane_stress_finite_deformation.dat";
  if (this->plane_strain)
    filename = "material_check_stress_plane_strain_finite_deformation.dat";
  else {
    SUCCEED();
    return;
  }

  this->initModel(_explicit_lumped_mass, filename);

  const auto & coordinates = this->mesh->getNodes();
  auto & displacement = this->model->getDisplacement();
  // set the position of all nodes to the static solution
  for (auto && tuple : zip(make_view(coordinates, this->dim),
                           make_view(displacement, this->dim))) {
    this->setLinearDOF(std::get<1>(tuple), std::get<0>(tuple));
  }
  for (UInt s = 0; s < 100; ++s) {
    this->model->solveStep();
  }

  auto ekin = this->model->getEnergy("kinetic");
  EXPECT_NEAR(0, ekin, 1e-16);

  this->checkAll();
}
