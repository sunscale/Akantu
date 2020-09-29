/**
 * @file   test_solid_mechanics_model_fixture.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Nov 14 2017
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Main solif mechanics test file
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
#include "communicator.hh"
#include "solid_mechanics_model.hh"
#include "test_gtest_utils.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_TEST_SOLID_MECHANICS_MODEL_FIXTURE_HH_
#define AKANTU_TEST_SOLID_MECHANICS_MODEL_FIXTURE_HH_

using namespace akantu;

// This fixture uses very small meshes with a volume of 1.
template <typename type_> class TestSMMFixture : public ::testing::Test {
public:
  static constexpr ElementType type = type_::value;
  static constexpr size_t spatial_dimension =
      ElementClass<type>::getSpatialDimension();

  void SetUp() override {
    this->mesh = std::make_unique<Mesh>(this->spatial_dimension);
    auto prank = Communicator::getStaticCommunicator().whoAmI();
    if (prank == 0) {
      this->mesh->read(this->mesh_file);
      if(spatial_dimension > 1 and mesh->getNbElement(spatial_dimension - 1) == 0) {
        MeshUtils::buildFacets(*this->mesh);
      }
    }

    mesh->distribute();

    SCOPED_TRACE(std::to_string(this->type).c_str());

    model = std::make_unique<SolidMechanicsModel>(*mesh, _all_dimensions,
                                                  std::to_string(this->type));
  }

  void initModel(const ID & input, const AnalysisMethod & analysis_method) {
    getStaticParser().parse(input);
    this->model->initFull(_analysis_method = analysis_method);

    if (analysis_method != _static) {
      auto time_step = this->model->getStableTimeStep() / 10.;
      this->model->setTimeStep(time_step);
    }
    // std::cout << "timestep: " << time_step << std::endl;

    if (this->dump_paraview) {
      std::stringstream base_name;
      base_name << "bar" << analysis_method << this->type;
      this->model->setBaseName(base_name.str());
      this->model->addDumpFieldVector("displacement");
      this->model->addDumpFieldVector("blocked_dofs");
      if (analysis_method != _static) {
        this->model->addDumpField("velocity");
        this->model->addDumpField("acceleration");
      }
      if (this->mesh->isDistributed()) {
        this->model->addDumpField("partitions");
      }
      this->model->addDumpFieldVector("external_force");
      this->model->addDumpFieldVector("internal_force");
      this->model->addDumpField("stress");
      this->model->addDumpField("strain");
    }
  }

  void TearDown() override {
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

protected:
  std::string mesh_file{std::to_string(this->type) + ".msh"};
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModel> model;
  bool dump_paraview{false};
};

template <typename type_> constexpr ElementType TestSMMFixture<type_>::type;
template <typename type_>
constexpr size_t TestSMMFixture<type_>::spatial_dimension;

template <typename T>
using is_not_pentahedron =
    aka::negation<aka::disjunction<is_element<T, _pentahedron_6>,
                                   is_element<T, _pentahedron_15>>>;

using TestElementTypesFiltered =
    tuple_filter_t<is_not_pentahedron, TestElementTypes>;

// using gtest_element_types = gtest_list_t<TestElementTypesFiltered>;
using gtest_element_types = gtest_list_t<TestElementTypes>;

TYPED_TEST_SUITE(TestSMMFixture, gtest_element_types);

#endif /* AKANTU_TEST_SOLID_MECHANICS_MODEL_FIXTURE_HH_ */
