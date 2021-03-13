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
#include "dof_manager.hh"
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

template <class Type>
class TestStructBernoulliDynamic : public TestStructuralFixture<Type> {
  using parent = TestStructuralFixture<Type>;

public:
  const UInt nb_element{40};
  const Real L{2};
  const Real le{L / nb_element};

  const UInt nb_nodes{nb_element + 1};
  const Real F{1e4};

  StructuralMaterial mat;

  void readMesh(std::string /*filename*/) override {
    MeshAccessor mesh_accessor(*this->mesh);
    auto & nodes = mesh_accessor.getNodes();
    nodes.resize(nb_nodes);

    nodes.set(0.);

    for (auto && data : enumerate(make_view(nodes, this->spatial_dimension))) {
      auto & node = std::get<1>(data);
      UInt i = std::get<0>(data);
      node[_x] = i * le;
    }

    this->mesh->addConnectivityType(parent::type);
    auto & connectivities = mesh_accessor.getConnectivity(parent::type);
    connectivities.resize(nb_element);
    for (auto && data : enumerate(make_view(connectivities, 2))) {
      UInt i = std::get<0>(data);
      auto & connectivity = std::get<1>(data);

      connectivity = {i, i + 1};
    }

    mesh_accessor.makeReady();
  }

  void setNormals() override {
    if (this->spatial_dimension != 3) {
      return;
    }
    auto & normals =
        this->mesh->template getData<Real>("extra_normal", parent::type);
    normals.resize(nb_element);

    for (auto && normal : make_view(normals, this->spatial_dimension)) {
      normal = {0., 0., 1.};
    }
  }

  AnalysisMethod getAnalysisMethod() const override {
    return _implicit_dynamic;
  }

  void addMaterials() override {
    this->mat.E = 1e9;
    this->mat.rho = 10;
    this->mat.I = 1;
    this->mat.Iz = 1;
    this->mat.Iy = 1;
    this->mat.A = 1;
    this->mat.GJ = 1;

    this->model->addMaterial(mat);
  }

  void setDirichletBCs() override {
    auto & boundary = this->model->getBlockedDOFs();
    boundary.set(false);

    boundary(0, _x) = true;
    boundary(0, _y) = true;

    boundary(nb_nodes - 1, _y) = true;

    if (this->spatial_dimension == 3) {
      boundary(0, _z) = true;
      boundary(nb_nodes - 1, _z) = true;
    }
  }

  void setNeumannBCs() override {
    auto node_to_print = nb_nodes / 2;
    // Forces
    auto & forces = this->model->getExternalForce();
    forces.zero();
    forces(node_to_print, _y) = F;
  }

  void assignMaterials() override {
    this->model->getElementMaterial(parent::type).set(0);
  }
};

using beam_types = gtest_list_t<std::tuple<element_type_t<_bernoulli_beam_2>,
                                           element_type_t<_bernoulli_beam_3>>>;

TYPED_TEST_SUITE(TestStructBernoulliDynamic, beam_types);

template <class Type>
void getElementMassMatrix(const StructuralMaterial & /*material*/, Real /*l*/,
                          Matrix<Real> & /*M*/) {}

template <class Type>
void getElementStifnessMatrix(const StructuralMaterial & /*material*/,
                              Real /*l*/, Matrix<Real> & /*M*/) {}

template <>
void getElementMassMatrix<element_type_t<_bernoulli_beam_2>>(
    const StructuralMaterial & material, Real l, Matrix<Real> & M) {
  auto A = material.A;
  auto rho = material.rho;
  // clang-format off
  M = rho * A * l / 420. * Matrix<Real>({
      {140,    0,      0,  70,     0,      0},
      {  0,  156,   22*l,   0,    54,  -13*l},
      {  0, 22*l,  4*l*l,   0,  13*l, -3*l*l},
      { 70,    0,      0, 140,     0,      0},
      {  0,   54,   13*l,   0,   156,  -22*l},
      {  0,-13*l, -3*l*l,   0, -22*l,  4*l*l}});
  // clang-format on
}

template <>
void getElementStifnessMatrix<element_type_t<_bernoulli_beam_2>>(
    const StructuralMaterial & material, Real l, Matrix<Real> & K) {
  auto E = material.E;
  auto A = material.A;
  auto I = material.I;

  auto l_2 = l * l;
  auto l_3 = l * l * l;
  // clang-format off
  K = Matrix<Real>({
      { E*A/l,           0,          0, -E*A/l,           0,          0},
      {     0,  12*E*I/l_3,  6*E*I/l_2,      0, -12*E*I/l_3,  6*E*I/l_2},
      {     0,   6*E*I/l_2,    4*E*I/l,      0,  -6*E*I/l_2,    2*E*I/l},
      {-E*A/l,           0,          0,  E*A/l,           0,          0},
      {     0, -12*E*I/l_3, -6*E*I/l_2,      0,  12*E*I/l_3, -6*E*I/l_2},
      {     0,   6*E*I/l_2,    2*E*I/l,      0,  -6*E*I/l_2,    4*E*I/l}});
  // clang-format on
}

template <>
void getElementMassMatrix<element_type_t<_bernoulli_beam_3>>(
    const StructuralMaterial & material, Real l, Matrix<Real> & M) {
  auto A = material.A;
  auto rho = material.rho;
  // clang-format off
  M = rho * A * l / 420. * Matrix<Real>({
      {140,     0,       0,   0,      0,      0,  70,     0,     0,   0,      0,      0},
      {  0,   156,       0,   0,      0,   22*l,   0,    54,     0,   0,      0,  -13*l},
      {  0,     0,     156,   0,  -22*l,      0,   0,     0,    54,   0,   13*l,      0},
      {  0,     0,       0, 140,      0,      0,   0,     0,     0,  70,      0,      0},
      {  0,     0,   -22*l,   0,  4*l*l,      0,   0,     0, -13*l,   0, -3*l*l,      0},
      {  0,  22*l,       0,   0,      0,  4*l*l,   0,  13*l,     0,   0,      0, -3*l*l},
      { 70,     0,       0,   0,      0,      0, 140,     0,     0,   0,      0,      0},
      {  0,    54,       0,   0,      0,   13*l,   0,   156,     0,   0,      0,  -22*l},
      {  0,     0,      54,   0,  -13*l,      0,   0,     0,   156,   0,   22*l,      0},
      {  0,     0,       0,  70,      0,      0,   0,     0,     0, 140,      0,      0},
      {  0,     0,    13*l,   0, -3*l*l,      0,   0,     0,  22*l,   0,  4*l*l,      0},
      {  0, -13*l,       0,   0,      0, -3*l*l,   0, -22*l,     0,   0,      0,  4*l*l}});
  // clang-format on
}

template <>
void getElementStifnessMatrix<element_type_t<_bernoulli_beam_3>>(
    const StructuralMaterial & material, Real l, Matrix<Real> & K) {
  auto E = material.E;
  auto A = material.A;
  auto Iy = material.Iy;
  auto Iz = material.Iz;
  auto GJ = material.GJ;

  auto a1 = E * A / l;
  auto b1 = 12 * E * Iz / l / l / l;
  auto b2 = 6 * E * Iz / l / l;
  auto b3 = 4 * E * Iz / l;
  auto b4 = 2 * E * Iz / l;

  auto c1 = 12 * E * Iy / l / l / l;
  auto c2 = 6 * E * Iy / l / l;
  auto c3 = 4 * E * Iy / l;
  auto c4 = 2 * E * Iy / l;

  auto d1 = GJ / l;

  // clang-format off
  K = Matrix<Real>({
      {  a1,   0,   0,   0,   0,   0, -a1,   0,   0,   0,   0,   0},
      {   0,  b1,   0,   0,   0,  b2,   0, -b1,   0,   0,   0,  b2},
      {   0,   0,  c1,   0, -c2,   0,   0,   0, -c1,   0, -c2,   0},
      {   0,   0,   0,  d1,   0,   0,   0,   0,   0, -d1,   0,   0},
      {   0,   0, -c2,   0,  c3,   0,   0,   0,  c2,   0,  c4,   0},
      {   0,  b2,   0,   0,   0,  b3,   0, -b2,   0,   0,   0,  b4},
      { -a1,   0,   0,   0,   0,   0,  a1,   0,   0,   0,   0,   0},
      {   0, -b1,   0,   0,   0, -b2,   0,  b1,   0,   0,   0, -b2},
      {   0,   0, -c1,   0,  c2,   0,   0,   0,  c1,   0,  c2,   0},
      {   0,   0,   0, -d1,   0,   0,   0,   0,   0,  d1,   0,   0},
      {   0,   0, -c2,   0,  c4,   0,   0,   0,  c2,   0,  c3,   0},
      {   0,  b2,   0,   0,   0,  b4,   0, -b2,   0,   0,   0,  b3}});
  // clang-format on
}

TYPED_TEST(TestStructBernoulliDynamic, TestBeamMatrices) {
  this->model->assembleMatrix("M");
  this->model->assembleMatrix("K");

  const auto & K = this->model->getDOFManager().getMatrix("K");
  const auto & M = this->model->getDOFManager().getMatrix("M");

  Matrix<Real> Ka(this->nb_nodes * this->ndof, this->nb_nodes * this->ndof, 0.);
  Matrix<Real> Ma(this->nb_nodes * this->ndof, this->nb_nodes * this->ndof, 0.);

  Matrix<Real> Ke(this->ndof * 2, this->ndof * 2);
  Matrix<Real> Me(this->ndof * 2, this->ndof * 2);

  getElementMassMatrix<TypeParam>(this->mat, this->le, Me);
  getElementStifnessMatrix<TypeParam>(this->mat, this->le, Ke);

  auto assemble = [&](auto && nodes, auto && M, auto && Me) {
    auto n1 = nodes[0];
    auto n2 = nodes[1];
    for (auto i : arange(this->ndof)) {
      for (auto j : arange(this->ndof)) {
        M(n1 * this->ndof + i, n1 * this->ndof + j) += Me(i, j);
        M(n2 * this->ndof + i, n2 * this->ndof + j) +=
            Me(this->ndof + i, this->ndof + j);
        M(n1 * this->ndof + i, n2 * this->ndof + j) += Me(i, this->ndof + j);
        M(n2 * this->ndof + i, n1 * this->ndof + j) += Me(this->ndof + i, j);
      }
    }
  };

  auto && connectivities = this->mesh->getConnectivity(this->type);

  for (auto && connectivity : make_view(connectivities, 2)) {
    assemble(connectivity, Ka, Ke);
    assemble(connectivity, Ma, Me);
  }

  auto tol = 1e-13;

  auto Ka_max = Ka.template norm<L_inf>();
  auto Ma_max = Ma.template norm<L_inf>();

  for (auto i : arange(Ka.rows())) {
    for (auto j : arange(Ka.cols())) {
      EXPECT_NEAR(Ka(i, j), K(i, j), tol * Ka_max);
      EXPECT_NEAR(Ma(i, j), M(i, j), tol * Ma_max);
    }
  }
}

TYPED_TEST(TestStructBernoulliDynamic, TestBeamOscilation) {
  Real time_step = 1e-6;
  this->model->setTimeStep(time_step);

  auto & solver = this->model->getNonLinearSolver();
  solver.set("max_iterations", 100);
  solver.set("threshold", 1e-8);
  solver.set("convergence_type", SolveConvergenceCriteria::_solution);

  auto node_to_print = this->nb_nodes / 2;
  auto & d = this->model->getDisplacement()(node_to_print, _y);

  std::ofstream pos;
  std::string filename = "position" + std::to_string(this->type) + ".csv";
  pos.open(filename);
  if (not pos.good()) {
    AKANTU_ERROR("Cannot open file \"position.csv\"");
  }
  pos << "id,time,position,solution" << std::endl;

//#define debug
#ifdef debug
  this->model->addDumpFieldVector("displacement");
  this->model->addDumpField("blocked_dofs");
  this->model->addDumpFieldVector("external_force");
  this->model->addDumpFieldVector("internal_force");
  this->model->addDumpFieldVector("acceleration");
  this->model->addDumpFieldVector("velocity");
  this->model->dump();
#endif

  this->model->getDisplacement().set(0.);

  Real tol = 1e-6;
  Real time = 0.;
  for (UInt s = 1; s < 300; ++s) {
    EXPECT_NO_THROW(this->model->solveStep());

    time = s * time_step;

    auto da = analytical_solution(time, this->L, this->mat.rho, this->mat.E,
                                  this->mat.A, this->mat.Iy, this->F);

    pos << s << "," << time << "," << d << "," << da << std::endl;

#ifdef debug
    this->model->dump();
#endif
    EXPECT_NEAR(d, da, tol);
  }

}
