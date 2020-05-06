/**
 * @file   patch_test_linear_anisotropic_explicit.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Till Junge <till.junge@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Tue Dec 05 2017
 * @date last modification: Tue Feb 13 2018
 *
 * @brief  patch test for elastic material in solid mechanics model
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
#include "non_linear_solver.hh"
#include "patch_test_linear_solid_mechanics_fixture.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

// Stiffness tensor, rotated by hand

/* -------------------------------------------------------------------------- */
TYPED_TEST(TestPatchTestSMMLinear, AnisotropicExplicit) {
  Real C[3][3][3][3] = {
      {{{112.93753505, 1.85842452538e-10, -4.47654358027e-10},
        {1.85847317471e-10, 54.2334345331, -3.69840984824},
        {-4.4764768395e-10, -3.69840984824, 56.848605217}},
       {{1.85847781609e-10, 25.429294233, -3.69840984816},
        {25.429294233, 3.31613847493e-10, -8.38797920011e-11},
        {-3.69840984816, -8.38804581349e-11, -1.97875715813e-10}},
       {{-4.47654358027e-10, -3.69840984816, 28.044464917},
        {-3.69840984816, 2.09374961813e-10, 9.4857455224e-12},
        {28.044464917, 9.48308098714e-12, -2.1367885239e-10}}},
      {{{1.85847781609e-10, 25.429294233, -3.69840984816},
        {25.429294233, 3.31613847493e-10, -8.38793479119e-11},
        {-3.69840984816, -8.38795699565e-11, -1.97876381947e-10}},
       {{54.2334345331, 3.31617400207e-10, 2.09372075233e-10},
        {3.3161562385e-10, 115.552705733, -3.15093728886e-10},
        {2.09372075233e-10, -3.15090176173e-10, 54.2334345333}},
       {{-3.69840984824, -8.38795699565e-11, 9.48219280872e-12},
        {-8.38795699565e-11, -3.1509195253e-10, 25.4292942335},
        {9.48441325477e-12, 25.4292942335, 3.69840984851}}},
      {{{-4.47653469848e-10, -3.69840984816, 28.044464917},
        {-3.69840984816, 2.09374073634e-10, 9.48752187924e-12},
        {28.044464917, 9.48552347779e-12, -2.1367885239e-10}},
       {{-3.69840984824, -8.3884899027e-11, 9.48219280872e-12},
        {-8.3884899027e-11, -3.150972816e-10, 25.4292942335},
        {9.48041645188e-12, 25.4292942335, 3.69840984851}},
       {{56.848605217, -1.97875493768e-10, -2.13681516925e-10},
        {-1.97877270125e-10, 54.2334345333, 3.69840984851},
        {-2.13683293282e-10, 3.69840984851, 112.93753505}}}};

  if (this->dim == 2) {
    for (UInt i = 0; i < this->dim; ++i) {
      for (UInt j = 0; j < this->dim; ++j) {
        for (UInt k = 0; k < this->dim; ++k) {
          for (UInt l = 0; l < this->dim; ++l) {
            C[i][j][k][l] = 0;
          }
        }
      }
    }
    C[0][0][0][0] = C[1][1][1][1] = 112.93753504999995;
    C[0][0][1][1] = C[1][1][0][0] = 51.618263849999984;
    C[0][1][0][1] = C[1][0][0][1] = C[0][1][1][0] = C[1][0][1][0] =
        22.814123549999987;
  }

  if (this->dim == 1) {
    C[0][0][0][0] = 105.092023;
  }
  this->initModel(_explicit_lumped_mass,
                  "material_anisotropic_" + std::to_string(this->dim) + ".dat");

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

  auto & mat = this->model->getMaterial(0);

  this->checkDOFs(displacement);
  this->checkGradient(mat.getGradU(this->type), displacement);

  this->result_tolerance = 1e-11;
  this->checkResults(
      [&](const Matrix<Real> & pstrain) {
        auto strain = (pstrain + pstrain.transpose()) / 2.;
        decltype(strain) stress(this->dim, this->dim);

        for (UInt i = 0; i < this->dim; ++i) {
          for (UInt j = 0; j < this->dim; ++j) {
            stress(i, j) = 0;
            for (UInt k = 0; k < this->dim; ++k) {
              for (UInt l = 0; l < this->dim; ++l) {
                stress(i, j) += C[i][j][k][l] * strain(k, l);
              }
            }
          }
        }
        return stress;
      },
      mat.getStress(this->type), displacement);
}

TYPED_TEST(TestPatchTestSMMLinear, AnisotropicStatic) {
  Real C[3][3][3][3] = {
      {{{112.93753505, 1.85842452538e-10, -4.47654358027e-10},
        {1.85847317471e-10, 54.2334345331, -3.69840984824},
        {-4.4764768395e-10, -3.69840984824, 56.848605217}},
       {{1.85847781609e-10, 25.429294233, -3.69840984816},
        {25.429294233, 3.31613847493e-10, -8.38797920011e-11},
        {-3.69840984816, -8.38804581349e-11, -1.97875715813e-10}},
       {{-4.47654358027e-10, -3.69840984816, 28.044464917},
        {-3.69840984816, 2.09374961813e-10, 9.4857455224e-12},
        {28.044464917, 9.48308098714e-12, -2.1367885239e-10}}},
      {{{1.85847781609e-10, 25.429294233, -3.69840984816},
        {25.429294233, 3.31613847493e-10, -8.38793479119e-11},
        {-3.69840984816, -8.38795699565e-11, -1.97876381947e-10}},
       {{54.2334345331, 3.31617400207e-10, 2.09372075233e-10},
        {3.3161562385e-10, 115.552705733, -3.15093728886e-10},
        {2.09372075233e-10, -3.15090176173e-10, 54.2334345333}},
       {{-3.69840984824, -8.38795699565e-11, 9.48219280872e-12},
        {-8.38795699565e-11, -3.1509195253e-10, 25.4292942335},
        {9.48441325477e-12, 25.4292942335, 3.69840984851}}},
      {{{-4.47653469848e-10, -3.69840984816, 28.044464917},
        {-3.69840984816, 2.09374073634e-10, 9.48752187924e-12},
        {28.044464917, 9.48552347779e-12, -2.1367885239e-10}},
       {{-3.69840984824, -8.3884899027e-11, 9.48219280872e-12},
        {-8.3884899027e-11, -3.150972816e-10, 25.4292942335},
        {9.48041645188e-12, 25.4292942335, 3.69840984851}},
       {{56.848605217, -1.97875493768e-10, -2.13681516925e-10},
        {-1.97877270125e-10, 54.2334345333, 3.69840984851},
        {-2.13683293282e-10, 3.69840984851, 112.93753505}}}};

  if (this->dim == 2) {
    for (UInt i = 0; i < this->dim; ++i) {
      for (UInt j = 0; j < this->dim; ++j) {
        for (UInt k = 0; k < this->dim; ++k) {
          for (UInt l = 0; l < this->dim; ++l) {
            C[i][j][k][l] = 0;
          }
        }
      }
    }
    C[0][0][0][0] = C[1][1][1][1] = 112.93753504999995;
    C[0][0][1][1] = C[1][1][0][0] = 51.618263849999984;
    C[0][1][0][1] = C[1][0][0][1] = C[0][1][1][0] = C[1][0][1][0] =
        22.814123549999987;
  }

  if (this->dim == 1) {
    C[0][0][0][0] = 105.092023;
  }

  this->initModel(_static,
                  "material_anisotropic_" + std::to_string(this->dim) + ".dat");

  auto & solver = this->model->getNonLinearSolver();
  solver.set("max_iterations", 2);
  solver.set("threshold", 2e-4);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);

  this->model->solveStep();

  auto & mat = this->model->getMaterial(0);

  const auto & displacement = this->model->getDisplacement();
  this->checkDOFs(displacement);
  this->checkGradient(mat.getGradU(this->type), displacement);

  this->result_tolerance = 1e-11;
  this->checkResults(
      [&](const Matrix<Real> & pstrain) {
        auto strain = (pstrain + pstrain.transpose()) / 2.;
        decltype(strain) stress(this->dim, this->dim);

        for (UInt i = 0; i < this->dim; ++i) {
          for (UInt j = 0; j < this->dim; ++j) {
            stress(i, j) = 0;
            for (UInt k = 0; k < this->dim; ++k) {
              for (UInt l = 0; l < this->dim; ++l) {
                stress(i, j) += C[i][j][k][l] * strain(k, l);
              }
            }
          }
        }
        return stress;
      },
      mat.getStress(this->type), displacement);
}
