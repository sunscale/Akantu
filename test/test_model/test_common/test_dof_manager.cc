/**
 * @file   test_dof_manager.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Wed Jan 30 2019
 *
 * @brief test the dof managers
 *
 * @section LICENSE
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
#include <dof_manager.hh>
#include <solver_vector_petsc.hh>
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <type_traits>
#include <numeric>
/* -------------------------------------------------------------------------- */

enum DOFManagerType { _dmt_default, _dmt_petsc };

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

/* -------------------------------------------------------------------------- */
using namespace akantu;
/* -------------------------------------------------------------------------- */
namespace akantu {
struct DOFManagerTester {
  DOFManagerTester(std::unique_ptr<DOFManager> dof_manager)
      : dof_manager(std::move(dof_manager)) {}

  DOFManager & operator*() { return *dof_manager; }
  DOFManager * operator->() { return dof_manager.get(); }
  void getArrayPerDOFs(const ID & id, SolverVector & vector,
                       Array<Real> & array) {
    dof_manager->getArrayPerDOFs(id, vector, array);
  }

  SolverVector & residual() { return *dof_manager->residual; }

  std::unique_ptr<DOFManager> dof_manager;
};
} // namespace akantu

template <class T> class DOFManagerFixture : public ::testing::Test {
public:
  constexpr static DOFManagerType type = T::value;

  void SetUp() override {
    mesh = std::make_unique<Mesh>(3);

    auto & communicator = Communicator::getStaticCommunicator();
    if (communicator.whoAmI() == 0)
      mesh->read("mesh.msh");
    mesh->distribute();

    nb_nodes = this->mesh->getNbNodes();
    nb_total_nodes = this->mesh->getNbGlobalNodes();

    auto && range_nodes = arange(nb_nodes);
    nb_pure_local = std::accumulate(range_nodes.begin(), range_nodes.end(), 0,
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

  decltype(auto) registerDOFs() {
    auto dof_manager = DOFManagerTester(this->alloc());

    this->dof1 = std::make_unique<Array<Real>>(nb_nodes, 3);
    dof_manager->registerDOFs("dofs1", *this->dof1, _dst_nodal);

    EXPECT_EQ(dof_manager.residual().size(), nb_total_nodes * 3);

    this->dof2 = std::make_unique<Array<Real>>(nb_pure_local, 5);
    dof_manager->registerDOFs("dofs2", *this->dof2, _dst_generic);

    EXPECT_EQ(dof_manager.residual().size(), nb_total_nodes * 8);
    return dof_manager;
  }

protected:
  Int nb_nodes{0}, nb_total_nodes{0}, nb_pure_local{0};
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<Array<Real>> dof1;
  std::unique_ptr<Array<Real>> dof2;
};

TYPED_TEST_CASE(DOFManagerFixture, dof_manager_types);
/* -------------------------------------------------------------------------- */

TYPED_TEST(DOFManagerFixture, Construction) {
  auto dof_manager = this->alloc();
}

TYPED_TEST(DOFManagerFixture, DoubleConstruction) {
  auto dof_manager = this->alloc();
  dof_manager = this->alloc();
}

TYPED_TEST(DOFManagerFixture, RegisterGenericDOF1) {
  auto dof_manager = this->alloc();

  Array<Real> dofs(this->nb_pure_local, 3);
  
  dof_manager->registerDOFs("dofs1", dofs, _dst_generic);
  EXPECT_GE(dof_manager.residual().size(), this->nb_total_nodes * 3);
}

TYPED_TEST(DOFManagerFixture, RegisterNodalDOF1) {
  auto dof_manager = this->alloc();

  Array<Real> dofs(this->nb_nodes, 3);
  dof_manager->registerDOFs("dofs1", dofs, _dst_nodal);
  EXPECT_GE(dof_manager.residual().size(), this->nb_total_nodes * 3);
}

TYPED_TEST(DOFManagerFixture, RegisterGenericDOF2) {
  auto dof_manager = this->alloc();

  Array<Real> dofs1(this->nb_pure_local, 3);
  dof_manager->registerDOFs("dofs1", dofs1, _dst_generic);
  EXPECT_EQ(dof_manager.residual().size(), this->nb_total_nodes * 3);

  Array<Real> dofs2(this->nb_pure_local, 5);
  dof_manager->registerDOFs("dofs2", dofs2, _dst_generic);
  EXPECT_EQ(dof_manager.residual().size(), this->nb_total_nodes * 8);
}

TYPED_TEST(DOFManagerFixture, RegisterNodalDOF2) {
  auto dof_manager = this->alloc();

  Array<Real> dofs1(this->nb_nodes, 3);
  dof_manager->registerDOFs("dofs1", dofs1, _dst_nodal);
  EXPECT_EQ(dof_manager.residual().size(), this->nb_total_nodes * 3);

  Array<Real> dofs2(this->nb_nodes, 5);
  dof_manager->registerDOFs("dofs2", dofs2, _dst_nodal);
  EXPECT_EQ(dof_manager.residual().size(), this->nb_total_nodes * 8);
}

TYPED_TEST(DOFManagerFixture, RegisterMixedDOF) {
  auto dof_manager = this->registerDOFs();
}

TYPED_TEST(DOFManagerFixture, Assemble) {
  auto dof_manager = this->registerDOFs();

  dof_manager.residual().clear();

  Array<Real> local1(*this->dof1);
  for(auto && data : enumerate(make_view(local1, local1.getNbComponent()))) {
    auto n = std::get<0>(data);
    auto & l = std::get<1>(data);
    l.set(1. * this->mesh->isLocalOrMasterNode(n));
  }

  dof_manager->assembleToResidual("dofs1", local1);

  Array<Real> local2(*this->dof2);
  local2.set(2.); 
  dof_manager->assembleToResidual("dofs2", local2);
  
  local1.set(0.);
  dof_manager.getArrayPerDOFs("dofs1", dof_manager.residual(), local1);
  for (auto && data : enumerate(make_view(local1, local1.getNbComponent()))) {
    if(this->mesh->isLocalOrMasterNode(std::get<0>(data))) {
      const auto & l = std::get<1>(data);
      auto e = (l - Vector<Real>{1., 1., 1.}).norm();
      ASSERT_EQ(e, 0.);
    }
  }
  
  local2.set(0.);
  dof_manager.getArrayPerDOFs("dofs2", dof_manager.residual(), local2);
  for (auto && l : make_view(local2, local2.getNbComponent())) {
    auto e = (l - Vector<Real>{2., 2., 2., 2., 2.}).norm();
    ASSERT_EQ(e, 0.);
  }
}
