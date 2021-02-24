/**
 * @file   test_structural_mechanics_model_fixture.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Tue Nov 14 2017
 * @date last modification: Fri Feb 09 2018
 *
 * @brief  Main test for structural model
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
#include "element_class_structural.hh"
#include "structural_mechanics_model.hh"
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_TEST_STRUCTURAL_MECHANICS_MODEL_FIXTURE_HH_
#define AKANTU_TEST_STRUCTURAL_MECHANICS_MODEL_FIXTURE_HH_

using namespace akantu;

template <typename type_> class TestStructuralFixture : public ::testing::Test {
public:
  static constexpr const ElementType type = type_::value;
  static constexpr const size_t spatial_dimension =
      ElementClass<type>::getSpatialDimension();
  static const UInt ndof = ElementClass<type>::getNbDegreeOfFreedom();

  void SetUp() override {
    const auto spatial_dimension = this->spatial_dimension;
    mesh = std::make_unique<Mesh>(spatial_dimension);
    readMesh(makeMeshName());

    std::stringstream element_type;
    element_type << this->type;
    model = std::make_unique<StructuralMechanicsModel>(*mesh, _all_dimensions,
                                                       element_type.str());

    setNormals();
    initModel();
  }

  virtual void initModel() {
    this->addMaterials();
    auto method = getAnalysisMethod();

    this->model->initFull(_analysis_method = method);

    this->assignMaterials();

    this->setDirichletBCs();
    this->setNeumannBCs();
  }

  virtual AnalysisMethod getAnalysisMethod() const { return _static; }

  virtual void readMesh(std::string filename) {
    mesh->read(filename, _miot_gmsh_struct);
  }

  virtual std::string makeMeshName() {
    std::stringstream element_type;
    element_type << type;
    SCOPED_TRACE(element_type.str().c_str());
    return element_type.str() + ".msh";
  }

  void TearDown() override {
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

  virtual void addMaterials() = 0;
  virtual void assignMaterials() = 0;
  virtual void setDirichletBCs() = 0;
  virtual void setNeumannBCs() = 0;
  virtual void setNormals() {}

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<StructuralMechanicsModel> model;
};

template <typename type_>
constexpr ElementType TestStructuralFixture<type_>::type;
template <typename type_>
constexpr size_t TestStructuralFixture<type_>::spatial_dimension;
template <typename type_> const UInt TestStructuralFixture<type_>::ndof;

using structural_types = gtest_list_t<StructuralTestElementTypes>;


#endif /* AKANTU_TEST_STRUCTURAL_MECHANICS_MODEL_FIXTURE_HH_ */
