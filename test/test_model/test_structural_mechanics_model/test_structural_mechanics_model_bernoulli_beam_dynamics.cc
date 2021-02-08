/**
 * @file   test_structural_mechanics_model_bernoulli_beam_dynamics.cc
 *
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 *
 * @date creation: Mon Jul 07 2014
 * @date last modification: Wed Feb 03 2016
 *
 * @brief  Test for _bernouilli_beam in dynamic
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
#include "test_structural_mechanics_model_fixture.hh"
/* -------------------------------------------------------------------------- */
#include "mesh_accessor.hh"
#include "non_linear_solver_newton_raphson.hh"
#include "structural_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>
#include <limits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
static Real analytical_solution(Real time, Real L, Real rho, Real E,
                                __attribute__((unused)) Real A, Real I,
                                Real F) {
  Real omega = M_PI * M_PI / L / L * sqrt(E * I / rho);
  Real sum = 0.;
  UInt i = 5;
  for (UInt n = 1; n <= i; n += 2) {
    sum += (1. - cos(n * n * omega * time)) / pow(n, 4);
  }

  return 2. * F * pow(L, 3) / pow(M_PI, 4) / E / I * sum;
}

class TestStructBernoulli3Dynamic
    : public TestStructuralFixture<element_type_t<_bernoulli_beam_3>> {
  using parent = TestStructuralFixture<element_type_t<_bernoulli_beam_3>>;

public:
  Real L{20};
  const UInt nb_element{10};
  UInt nb_nodes;
  StructuralMaterial mat;
  const Real F{1.};

  void readMesh(std::string /*filename*/) override {
    nb_nodes = nb_element + 1;
    auto length = L / nb_element;
    MeshAccessor mesh_accessor(*this->mesh);
    auto & nodes = mesh_accessor.getNodes();
    nodes.resize(nb_nodes);

    this->mesh->addConnectivityType(_bernoulli_beam_3);
    auto & connectivities = mesh_accessor.getConnectivity(parent::type);
    connectivities.resize(nb_element);

    this->mesh->getElementalData<Real>("extra_normal")
        .initialize(*this->mesh, _element_kind = _ek_structural,
                    _nb_component = 3, _with_nb_element = true,
                    _default_value = 0.);

    auto & normals = this->mesh->getData<Real>("extra_normal", parent::type);
    normals.resize(nb_element);

    for (auto && data : enumerate(make_view(nodes, 3))) {
      auto & node = std::get<1>(data);
      UInt i = std::get<0>(data);
      node = {i * length, 0., 0.};
    }

    for (auto && data :
         enumerate(make_view(connectivities, 2), make_view(normals, 3))) {
      UInt i = std::get<0>(data);
      auto & connectivity = std::get<1>(data);
      auto & normal = std::get<2>(data);

      connectivity = {i, i + 1};
      normal = {0., 0., 1.};
    }

    mesh_accessor.makeReady();
  }

  AnalysisMethod getAnalysisMethod() const override {
    return _implicit_dynamic;
  }

  void addMaterials() override {
    this->mat.E = 1e9;
    this->mat.rho = 1;
    this->mat.Iz = 1;
    this->mat.Iy = 1;
    this->mat.A = 1;
    this->mat.GJ = 1;

    this->model->addMaterial(mat);
  }

  void setDirichlets() override {
    auto boundary = this->model->getBlockedDOFs().begin(parent::ndof);
    // clang-format off
    Vector<bool> boundary_b = boundary[0];
    Vector<bool> boundary_e = boundary[nb_nodes];
    boundary_b = {true, true, true, false, false, false};
    boundary_e = {false, true, true, false, false, false};
    // clang-format on
  }

  void setNeumanns() override {
    auto node_to_print = nb_nodes / 2 + 1;
    // Forces
    auto & forces = this->model->getExternalForce();
    forces(node_to_print - 1, _y) = F;
  }

  void assignMaterials() override {
    model->getElementMaterial(parent::type).set(0);
  }
};

TEST_F(TestStructBernoulli3Dynamic, TestBeamOscilation) {
  auto & solver = this->model->getNonLinearSolver();
  solver.set("max_iterations", 100);
  solver.set("threshold", 1e-8);
  solver.set("convergence_type", SolveConvergenceCriteria::_solution);

  auto node_to_print = this->nb_nodes / 2 + 1;
  auto & d = this->model->getDisplacement()(node_to_print, 1);

  Real time_step = 1e-5;
  this->model->setTimeStep(time_step);


  std::ofstream pos;
  pos.open("position.csv");
  if (not pos.good()) {
    AKANTU_ERROR("Cannot open file \"position.csv\"");
  }
  pos << "id,time,position,solution" << std::endl;

  Real tol = 1e-10;
  Real time = 0.;
  for (UInt s = 1; time < 0.1606; ++s) {
    EXPECT_NO_THROW(this->model->solveStep());

    time = s * time_step;

    auto da = analytical_solution(time, this->L, this->mat.rho, this->mat.E,
                                  this->mat.A, this->mat.Iy, this->F);

    pos << s << "," << time << "," << d << "," << da <<std::endl;
    //EXPECT_NEAR(d , da, tol);
  }
}
