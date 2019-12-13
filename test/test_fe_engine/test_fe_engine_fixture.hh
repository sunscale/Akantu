/**
 * @file   test_fe_engine_fixture.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Nov 14 2017
 * @date last modification: Mon Feb 19 2018
 *
 * @brief  Fixture for feengine tests
 *
 * @section LICENSE
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
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <element_class.hh>
#include <fe_engine.hh>
#include <integrator_gauss.hh>
#include <shape_lagrange.hh>
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TEST_FE_ENGINE_FIXTURE_HH__
#define __AKANTU_TEST_FE_ENGINE_FIXTURE_HH__

using namespace akantu;

/// Generic class for FEEngine tests
template <typename type_, template <ElementKind> class shape_t,
          ElementKind kind = _ek_regular>
class TestFEMBaseFixture : public ::testing::Test {
public:
  static constexpr const ElementType type = type_::value;
  static constexpr const size_t dim = ElementClass<type>::getSpatialDimension();
  using FEM = FEEngineTemplate<IntegratorGauss, shape_t, kind>;

  /// Setup reads mesh corresponding to element type and initializes an FEEngine
  void SetUp() override {
    const auto dim = this->dim;
    mesh = std::make_unique<Mesh>(dim);

    std::stringstream meshfilename;
    meshfilename << type << ".msh";
    this->readMesh(meshfilename.str());

    lower = mesh->getLowerBounds();
    upper = mesh->getUpperBounds();

    nb_element = this->mesh->getNbElement(type);

    fem = std::make_unique<FEM>(*mesh, dim, "my_fem");
    nb_quadrature_points_total =
        GaussIntegrationElement<type>::getNbQuadraturePoints() * nb_element;

    SCOPED_TRACE(std::to_string(type));
  }

  void TearDown() override {
    fem.reset(nullptr);
    mesh.reset(nullptr);
  }

  /// Should be reimplemented if further treatment of the mesh is needed
  virtual void readMesh(std::string file_name) { mesh->read(file_name); }

protected:
  std::unique_ptr<FEM> fem;
  std::unique_ptr<Mesh> mesh;
  UInt nb_element;
  UInt nb_quadrature_points_total;
  Vector<Real> lower;
  Vector<Real> upper;
};

template <typename type_, template <ElementKind> class shape_t,
          ElementKind kind>
constexpr const ElementType TestFEMBaseFixture<type_, shape_t, kind>::type;

template <typename type_, template <ElementKind> class shape_t,
          ElementKind kind>
constexpr const size_t TestFEMBaseFixture<type_, shape_t, kind>::dim;

/* -------------------------------------------------------------------------- */
/// Base class for test with Lagrange FEEngine and regular elements
template <typename type_>
using TestFEMFixture = TestFEMBaseFixture<type_, ShapeLagrange, _ek_regular>;

/* -------------------------------------------------------------------------- */

using fe_engine_types = gtest_list_t<TestElementTypes>;

TYPED_TEST_SUITE(TestFEMFixture, fe_engine_types);

#endif /* __AKANTU_TEST_FE_ENGINE_FIXTURE_HH__ */
