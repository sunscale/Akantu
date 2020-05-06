/**
 * @file   patch_test_linear_fixture.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jan 30 2018
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Fixture for linear patch tests
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
#include "element_group.hh"
#include "mesh_utils.hh"
#include "model.hh"
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PATCH_TEST_LINEAR_FIXTURE_HH__
#define __AKANTU_PATCH_TEST_LINEAR_FIXTURE_HH__

//#define DEBUG_TEST

using namespace akantu;

template <typename type_, typename M>
class TestPatchTestLinear : public ::testing::Test {
public:
  static constexpr ElementType type = type_::value;
  static constexpr size_t dim = ElementClass<type>::getSpatialDimension();

  virtual void SetUp() {
    mesh = std::make_unique<Mesh>(dim);
    mesh->read(std::to_string(type) + ".msh");
    MeshUtils::buildFacets(*mesh);
    mesh->createBoundaryGroupFromGeometry();

    model = std::make_unique<M>(*mesh);
  }

  virtual void TearDown() {
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

  virtual void initModel(const AnalysisMethod & method,
                         const std::string & material_file) {
    debug::setDebugLevel(dblError);
    getStaticParser().parse(material_file);

    this->model->initFull(_analysis_method = method);
    this->applyBC();

    if (method != _static)
      this->model->setTimeStep(0.8 * this->model->getStableTimeStep());
  }

  virtual void applyBC() {
    auto & boundary = this->model->getBlockedDOFs();

    for (auto & eg : mesh->iterateElementGroups()) {
      for (const auto & node : eg.getNodeGroup()) {
        for (UInt s = 0; s < boundary.getNbComponent(); ++s) {
          boundary(node, s) = true;
        }
      }
    }
  }

  virtual void applyBConDOFs(const Array<Real> & dofs) {
    const auto & coordinates = this->mesh->getNodes();
    for (auto & eg : this->mesh->iterateElementGroups()) {
      for (const auto & node : eg.getNodeGroup()) {
        this->setLinearDOF(dofs.begin(dofs.getNbComponent())[node],
                           coordinates.begin(this->dim)[node]);
      }
    }
  }

  template <typename V> Matrix<Real> prescribed_gradient(const V & dof) {
    Matrix<Real> gradient(dof.getNbComponent(), dim);

    for (UInt i = 0; i < gradient.rows(); ++i) {
      for (UInt j = 0; j < gradient.cols(); ++j) {
        gradient(i, j) = alpha(i, j + 1);
      }
    }
    return gradient;
  }

  template <typename Gradient, typename DOFs>
  void checkGradient(const Gradient & gradient, const DOFs & dofs) {
    auto pgrad = prescribed_gradient(dofs);

    for (auto & grad :
         make_view(gradient, gradient.getNbComponent() / dim, dim)) {
      auto diff = grad - pgrad;
      auto gradient_error =
          diff.template norm<L_inf>() / grad.template norm<L_inf>();

      EXPECT_NEAR(0, gradient_error, gradient_tolerance);
    }
  }

  template <typename presult_func_t, typename Result, typename DOFs>
  void checkResults(presult_func_t && presult_func, const Result & results,
                    const DOFs & dofs) {
    auto presult = presult_func(prescribed_gradient(dofs));
    for (auto & result :
         make_view(results, results.getNbComponent() / dim, dim)) {
      auto diff = result - presult;
      auto result_error =
          diff.template norm<L_inf>() / presult.template norm<L_inf>();

      EXPECT_NEAR(0, result_error, result_tolerance);
    }
  }

  template <typename V1, typename V2>
  void setLinearDOF(V1 && dof, V2 && coord) {
    for (UInt i = 0; i < dof.size(); ++i) {
      dof(i) = this->alpha(i, 0);
      for (UInt j = 0; j < coord.size(); ++j) {
        dof(i) += this->alpha(i, j + 1) * coord(j);
      }
    }
  }

  template <typename V> void checkDOFs(V && dofs) {
    const auto & coordinates = mesh->getNodes();
    Vector<Real> ref_dof(dofs.getNbComponent());

    for (auto && tuple : zip(make_view(coordinates, dim),
                             make_view(dofs, dofs.getNbComponent()))) {
      setLinearDOF(ref_dof, std::get<0>(tuple));
      auto diff = std::get<1>(tuple) - ref_dof;
      auto dofs_error = diff.template norm<L_inf>();

      EXPECT_NEAR(0, dofs_error, dofs_tolerance);
    }
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<M> model;
  Matrix<Real> alpha{{0.01, 0.02, 0.03, 0.04},
                     {0.05, 0.06, 0.07, 0.08},
                     {0.09, 0.10, 0.11, 0.12}};

  Real gradient_tolerance{1e-13};
  Real result_tolerance{1e-13};
  Real dofs_tolerance{1e-15};
};

template <typename type_, typename M>
constexpr ElementType TestPatchTestLinear<type_, M>::type;

template <typename tuple_, typename M>
constexpr size_t TestPatchTestLinear<tuple_, M>::dim;

#endif /* __AKANTU_PATCH_TEST_LINEAR_FIXTURE_HH__ */
