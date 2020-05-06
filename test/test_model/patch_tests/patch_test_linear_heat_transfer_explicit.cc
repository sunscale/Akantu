/**
 * @file   patch_test_linear_heat_transfer_explicit.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jan 30 2018
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  HeatTransfer patch test
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
#include "patch_test_linear_heat_transfer_fixture.hh"
/* -------------------------------------------------------------------------- */

TYPED_TEST(TestPatchTestHTMLinear, Explicit) {
  this->initModel(_explicit_lumped_mass, "heat_transfer_input.dat");

  const auto & coordinates = this->mesh->getNodes();
  auto & temperature = this->model->getTemperature();
  // set the position of all nodes to the static solution
  for (auto && tuple :
       zip(make_view(coordinates, this->dim), make_view(temperature, 1))) {
    this->setLinearDOF(std::get<1>(tuple), std::get<0>(tuple));
  }

  for (UInt s = 0; s < 100; ++s) {
    this->model->solveStep();
  }

  this->checkAll();
}
