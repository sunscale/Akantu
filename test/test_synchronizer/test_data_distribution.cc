/**
 * @file   test_data_distribution.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Sep 05 2014
 * @date last modification: Fri Jan 26 2018
 *
 * @brief  Test the mesh distribution on creation of a distributed synchonizer
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "test_synchronizers_fixture.hh"
/* -------------------------------------------------------------------------- */

TEST_F(TestSynchronizerFixture, DataDistribution) {
  auto & barycenters = this->mesh->getElementalData<Real>("barycenters");
  auto spatial_dimension = this->mesh->getSpatialDimension();
  barycenters.initialize(*this->mesh, _spatial_dimension = _all_dimensions,
                         _nb_component = spatial_dimension);

  this->initBarycenters(barycenters, *this->mesh);

  auto & gids = this->mesh->getNodalData<UInt>("gid");
  gids.resize(this->mesh->getNbNodes());

  for (auto && data : enumerate(gids)) {
    std::get<1>(data) = std::get<0>(data);
  }

  this->distribute();

  for (auto && ghost_type : ghost_types) {
    for (const auto & type :
         this->mesh->elementTypes(_all_dimensions, ghost_type)) {
      auto & barycenters =
          this->mesh->getData<Real>("barycenters", type, ghost_type);

      for (auto && data :
           enumerate(make_view(barycenters, spatial_dimension))) {
        Element element{type, UInt(std::get<0>(data)), ghost_type};
        Vector<Real> barycenter(spatial_dimension);
        this->mesh->getBarycenter(element, barycenter);

        auto dist = (std::get<1>(data) - barycenter).template norm<L_inf>();
        EXPECT_NEAR(dist, 0, 1e-7);
      }
    }
  }

  if (psize > 1) {
    for (auto && data : zip(gids, this->mesh->getGlobalNodesIds())) {
      EXPECT_EQ(std::get<0>(data), std::get<1>(data));
    }
  }
}

TEST_F(TestSynchronizerFixture, DataDistributionTags) {
  this->distribute();

  for (const auto & type : this->mesh->elementTypes(_all_dimensions)) {
    auto & tags = this->mesh->getData<UInt>("tag_0", type);
    Array<UInt>::const_vector_iterator tags_it = tags.begin(1);
    Array<UInt>::const_vector_iterator tags_end = tags.end(1);

    // The number of tags should match the number of elements on rank"
    EXPECT_EQ(this->mesh->getNbElement(type), tags.size());
  }
}
