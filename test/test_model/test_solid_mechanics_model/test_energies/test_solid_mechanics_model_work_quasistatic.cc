/**
 * @file   test_solid_mechanics_model_work_quasistatic.cc
 *
 * @author Tobias Brink <tobias.brink@epfl.ch>
 *
 * @date creation: Wed Nov 29 2017
 * @date last modification: Fri Jan 26 2018
 *
 * @brief  test work in quasistatic
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
 * @section description
 *
 * Assuming that the potential energy of a linear elastic material
 * works correctly, the work in a static simulation must equal the
 * potential energy of the material. Since the work in static is an
 * infinitesimal work Fds, we need to integrate by increasing F from 0
 * to F_final in steps. This test uses one Dirichlet boundary
 * condition (with u = 0.0, 0.1, and -0.1) and one Neumann boundary
 * condition for F on the opposite side. The final work must be the
 * same for all u.
 *
 */

/* -------------------------------------------------------------------------- */
#include "../test_solid_mechanics_model_fixture.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

namespace {

TYPED_TEST(TestSMMFixture, WorkQuasistatic) {
  const auto spatial_dimension = this->spatial_dimension;

  getStaticParser().parse("test_solid_mechanics_model_"
                          "work_material.dat");

  /// model initialization
  this->model->initFull(_analysis_method = _static);

  /// Create a node group for Neumann BCs.
  auto & apply_force_grp = this->mesh->createNodeGroup("apply_force");
  auto & fixed_grp = this->mesh->createNodeGroup("fixed");

  const auto & pos = this->mesh->getNodes();
  auto & flags = this->model->getBlockedDOFs();
  auto & lower = this->mesh->getLowerBounds();
  auto & upper = this->mesh->getUpperBounds();

  UInt i = 0;
  for (auto && data : zip(make_view(pos, spatial_dimension),
                          make_view(flags, spatial_dimension))) {
    const auto & posv = std::get<0>(data);
    auto & flag = std::get<1>(data);
    if (posv(_x) > upper(_x) - 1e-6) {
      apply_force_grp.add(i);
    } else if (posv(_x) < lower(_x) + 1e-6) {
      fixed_grp.add(i);

      if ((spatial_dimension > 1) and (posv(_y) < lower(_y) + 1e-6)) {
        flag(_y) = true;
        if ((spatial_dimension > 2) and (posv(_z) < lower(_z) + 1e-6)) {
          flag(_z) = true;
        }
      }
    }
    ++i;
  }

  this->mesh->createElementGroupFromNodeGroup("el_apply_force", "apply_force",
                                              spatial_dimension - 1);
  this->mesh->createElementGroupFromNodeGroup("el_fixed", "fixed",
                                              spatial_dimension - 1);

  std::vector<Real> displacements{0.0, 0.1, -0.1};
  for (auto && u : displacements) {
    this->model->applyBC(BC::Dirichlet::FixedValue(u, _x), "el_fixed");

    Vector<Real> surface_traction(spatial_dimension);
    Real work = 0.0;
    Real Epot;
    static const UInt N = 100;
    for (UInt i = 0; i <= N; ++i) {
      this->model->getExternalForce().zero(); // reset external forces to zero

      surface_traction(_x) = (1.0 * i) / N;

      if (spatial_dimension == 1) {
        // \TODO: this is a hack to work
        //      around non-implemented
        //      BC::Neumann::FromTraction for 1D
        auto & force = this->model->getExternalForce();
        for (auto && pair : zip(make_view(pos, spatial_dimension),
                                make_view(force, spatial_dimension))) {
          auto & posv = std::get<0>(pair);
          auto & forcev = std::get<1>(pair);
          if (posv(_x) > upper(_x) - 1e-6) {
            forcev(_x) = surface_traction(_x);
          }
        }
      } else {
        this->model->applyBC(BC::Neumann::FromTraction(surface_traction),
                             "el_apply_force");
      }

      /// Solve.
      this->model->solveStep();

      Epot = this->model->getEnergy("potential");
      // In static, this is infinitesimal work!
      auto Fds = this->model->getEnergy("external work");

      work += Fds; // integrate

      /// Check that no work was done for zero force.
      if (i == 0) {
        EXPECT_NEAR(work, 0.0, 1e-12);
      }
    }

    // Due to the finite integration steps, we make a rather large error
    // in our work integration, thus the allowed delta is 1e-2.
    EXPECT_NEAR(work, Epot, 1e-2);
  }
}

} // namespace
