/**
 * @file   test_synchronizers_fixture.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jan 26 2018
 *
 * @brief  Fixture for synchronizer tests
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
#include "aka_iterators.hh"
#include "communicator.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

using namespace akantu;

class TestSynchronizerFixture : public ::testing::Test {
public:
  virtual void SetUp() {
    const UInt spatial_dimension = 3;

    mesh = std::make_unique<Mesh>(spatial_dimension);

    const auto & comm = Communicator::getStaticCommunicator();
    prank = comm.whoAmI();
    psize = comm.getNbProc();

    if (prank == 0) {
      this->mesh->read("cube.msh");
    }
  }

  virtual void TearDown() { this->mesh.reset(nullptr); }

  void initBarycenters(ElementTypeMapArray<Real> & barycenters, Mesh & mesh) {
    auto spatial_dimension = mesh.getSpatialDimension();
    barycenters.initialize(mesh, _spatial_dimension = _all_dimensions,
                           _nb_component = spatial_dimension,
                           _with_nb_element = true);

    for (auto && ghost_type : ghost_types) {
      for (const auto & type : mesh.elementTypes(_all_dimensions, ghost_type)) {
        for (auto && data : enumerate(
                 make_view(barycenters(type, ghost_type), spatial_dimension))) {
          Element element{type, UInt(std::get<0>(data)), ghost_type};
          mesh.getBarycenter(element, std::get<1>(data));
        }
      }
    }
  }

  void distribute() { this->mesh->distribute(); }

protected:
  std::unique_ptr<Mesh> mesh;
  Int prank;
  Int psize;
};
