/**
 * @file   test_dof_manager.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Wed Jan 30 2019
 *
 * @brief test the dof managers
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <dof_manager.hh>
#include <mesh_partition_scotch.hh>
#include <mesh_utils.hh>
/* -------------------------------------------------------------------------- */

#include <gtest/gtest.h>
#include <numeric>
#include <string>
#include <type_traits>
/* -------------------------------------------------------------------------- */
namespace akantu {
enum DOFManagerType { _dmt_default, _dmt_petsc };
}
AKANTU_ENUM_HASH(DOFManagerType)

using namespace akantu;

// defined as struct to get there names in gtest outputs
struct _dof_manager_default
    : public std::integral_constant<DOFManagerType, _dmt_default> {};
struct _dof_manager_petsc
    : public std::integral_constant<DOFManagerType, _dmt_petsc> {};

using dof_manager_types = ::testing::Types<
#ifdef AKANTU_USE_PETSC
    _dof_manager_petsc,
#endif
    _dof_manager_default>;

namespace std {

std::string to_string(const DOFManagerType & type) {
  std::unordered_map<DOFManagerType, std::string> map{
#ifdef AKANTU_USE_PETSC
      {_dmt_petsc, "petsc"},
#endif
      {_dmt_default, "default"},
  };
  return map.at(type);
}

} // namespace std

/* -------------------------------------------------------------------------- */
using namespace akantu;
/* -------------------------------------------------------------------------- */
namespace akantu {
class DOFManagerTester {
public:
  DOFManagerTester(std::unique_ptr<DOFManager> dof_manager)
      : dof_manager(std::move(dof_manager)) {}

  DOFManager & operator*() { return *dof_manager; }
  DOFManager * operator->() { return dof_manager.get(); }
  void getArrayPerDOFs(const ID & id, SolverVector & vector,
                       Array<Real> & array) {
    dof_manager->getArrayPerDOFs(id, vector, array);
  }

  SolverVector & residual() { return *dof_manager->residual; }

private:
  std::unique_ptr<DOFManager> dof_manager;
};
} // namespace akantu

template <class T> class DOFManagerFixture : public ::testing::Test {
public:
  constexpr static DOFManagerType type = T::value;
  constexpr static UInt dim = 3;
  void SetUp() override {
    mesh = std::make_unique<Mesh>(this->dim);

    auto & communicator = Communicator::getStaticCommunicator();
    if (communicator.whoAmI() == 0) {
      mesh->read("mesh.msh");
    }
    mesh->distribute();

    nb_nodes = this->mesh->getNbNodes();
    nb_total_nodes = this->mesh->getNbGlobalNodes();

    auto && range_nodes = arange(nb_nodes);
    nb_pure_local =
        std::accumulate(range_nodes.begin(), range_nodes.end(), 0,
                        [&](auto && init, auto && val) {
                          return init + mesh->isLocalOrMasterNode(val);
                        });
  }
  void TearDown() override {
    mesh.reset();
    dof1.reset();
    dof2.reset();
  }

  decltype(auto) alloc() {
    std::unordered_map<DOFManagerType, std::string> types{
      {_dmt_default, "default"}, {_dmt_petsc, "petsc"}};

    return DOFManagerTester(DOFManagerFactory::getInstance().allocate(
        types[T::value], *mesh, "dof_manager", 0));
  }

  decltype(auto) registerDOFs(DOFSupportType dst1, DOFSupportType dst2) {
    auto dof_manager = DOFManagerTester(this->alloc());

    auto n1 = dst1 == _dst_nodal ? nb_nodes : nb_pure_local;
    this->dof1 = std::make_unique<Array<Real>>(n1, 3);

    dof_manager->registerDOFs("dofs1", *this->dof1, dst1);

    EXPECT_EQ(dof_manager.residual().size(), nb_total_nodes * 3);

    auto n2 = dst2 == _dst_nodal ? nb_nodes : nb_pure_local;
    this->dof2 = std::make_unique<Array<Real>>(n2, 5);

    dof_manager->registerDOFs("dofs2", *this->dof2, dst2);

    EXPECT_EQ(dof_manager.residual().size(), nb_total_nodes * 8);
    return dof_manager;
  }

protected:
  Int nb_nodes{0}, nb_total_nodes{0}, nb_pure_local{0};
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<Array<Real>> dof1;
  std::unique_ptr<Array<Real>> dof2;
};

template <class T> constexpr DOFManagerType DOFManagerFixture<T>::type;
template <class T> constexpr UInt DOFManagerFixture<T>::dim;

TYPED_TEST_SUITE(DOFManagerFixture, dof_manager_types);

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, Construction) {
  auto dof_manager = this->alloc();
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, DoubleConstruction) {
  auto dof_manager = this->alloc();
  dof_manager = this->alloc();
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, RegisterGenericDOF1) {
  auto dof_manager = this->alloc();

  Array<Real> dofs(this->nb_pure_local, 3);

  dof_manager->registerDOFs("dofs1", dofs, _dst_generic);
  EXPECT_GE(dof_manager.residual().size(), this->nb_total_nodes * 3);
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, RegisterNodalDOF1) {
  auto dof_manager = this->alloc();

  Array<Real> dofs(this->nb_nodes, 3);
  dof_manager->registerDOFs("dofs1", dofs, _dst_nodal);
  EXPECT_GE(dof_manager.residual().size(), this->nb_total_nodes * 3);
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, RegisterGenericDOF2) {
  this->registerDOFs(_dst_generic, _dst_generic);
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, RegisterNodalDOF2) {
  this->registerDOFs(_dst_nodal, _dst_nodal);
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, RegisterMixedDOF) {
  auto dof_manager = this->registerDOFs(_dst_nodal, _dst_generic);
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, AssembleVector) {
  auto dof_manager = this->registerDOFs(_dst_nodal, _dst_generic);

  dof_manager.residual().clear();

  for (auto && data :
       enumerate(make_view(*this->dof1, this->dof1->getNbComponent()))) {
    auto n = std::get<0>(data);
    auto & l = std::get<1>(data);
    l.set(1. * this->mesh->isLocalOrMasterNode(n));
  }

  this->dof2->set(2.);

  dof_manager->assembleToResidual("dofs1", *this->dof1);
  dof_manager->assembleToResidual("dofs2", *this->dof2);

  this->dof1->set(0.);
  this->dof2->set(0.);

  dof_manager.getArrayPerDOFs("dofs1", dof_manager.residual(), *this->dof1);
  for (auto && data :
       enumerate(make_view(*this->dof1, this->dof1->getNbComponent()))) {
    if (this->mesh->isLocalOrMasterNode(std::get<0>(data))) {
      const auto & l = std::get<1>(data);
      auto e = (l - Vector<Real>{1., 1., 1.}).norm();
      ASSERT_EQ(e, 0.);
    }
  }

  dof_manager.getArrayPerDOFs("dofs2", dof_manager.residual(), *this->dof2);
  for (auto && l : make_view(*this->dof2, this->dof2->getNbComponent())) {
    auto e = (l - Vector<Real>{2., 2., 2., 2., 2.}).norm();
    ASSERT_EQ(e, 0.);
  }
}

/* -------------------------------------------------------------------------- */
TYPED_TEST(DOFManagerFixture, AssembleMatrixNodal) {
  auto dof_manager = this->registerDOFs(_dst_nodal, _dst_nodal);

  auto && K = dof_manager->getNewMatrix("K", _symmetric);
  K.clear();

  auto && elemental_matrix = std::make_unique<Array<Real>>(
      this->mesh->getNbElement(this->dim), 8 * 3 * 8 * 3);

  for (auto && m : make_view(*elemental_matrix, 8 * 3, 8 * 3)) {
    m.set(1.);
  }

  dof_manager->assembleElementalMatricesToMatrix(
      "K", "dofs1", *elemental_matrix, _hexahedron_8);

  elemental_matrix = std::make_unique<Array<Real>>(
      this->mesh->getNbElement(this->dim), 8 * 5 * 8 * 5);

  for (auto && m : make_view(*elemental_matrix, 8 * 5, 8 * 5)) {
    m.set(1.);
  }

  dof_manager->assembleElementalMatricesToMatrix(
      "K", "dofs2", *elemental_matrix, _hexahedron_8);

  CSR<Element> node_to_elem;
  MeshUtils::buildNode2Elements(*this->mesh, node_to_elem, this->dim);

  dof_manager.residual().clear();

  for (auto && data :
       enumerate(zip(make_view(*this->dof1, this->dof1->getNbComponent()),
                     make_view(*this->dof2, this->dof2->getNbComponent())))) {
    auto n = std::get<0>(data);
    auto & l1 = std::get<0>(std::get<1>(data));
    auto & l2 = std::get<1>(std::get<1>(data));
    auto v = 1. * this->mesh->isLocalOrMasterNode(n);
    l1.set(v);
    l2.set(v);
  }

  dof_manager->assembleToResidual("dofs1", *this->dof1);
  dof_manager->assembleToResidual("dofs2", *this->dof2);

  for (auto && n : arange(this->nb_nodes)) {
    if (not this->mesh->isLocalOrMasterNode(n)) {
    }
  }
}
