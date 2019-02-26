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
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

enum DOFManagerType { _dmt_default, _dmt_petsc };

// defined as struct to get there names in gtest outputs
struct _dof_manager_default
    : public std::integral_constant<DOFManagerType, _dmt_default> {};
struct _dof_manager_petsc
    : public std::integral_constant<DOFManagerType, _dmt_petsc> {};

using dof_manager_types = ::testing::Types<
    _dof_manager_default
#ifdef AKANTU_USE_PETSC
    , _dof_manager_petsc
#endif
  >;

/* -------------------------------------------------------------------------- */
using namespace akantu;
/* -------------------------------------------------------------------------- */

template <class T> class DOFManagerFixture : public ::testing::Test {
public:
  void SetUp() override {
    mesh = std::make_unique<Mesh>(3);

    auto & communicator = Communicator::getStaticCommunicator();
    if (communicator.whoAmI() == 0)
      mesh->read("mesh.msh");
    mesh->distribute();
    
    dof1 = std::make_unique<Array<Real>>(mesh->getNbNodes(), 3);
    dof2 = std::make_unique<Array<Real>>(mesh->getNbNodes(), 5);
  }
  void TearDown() override {
    mesh.reset();
    dof1.reset();
    dof2.reset();
  }

  decltype(auto) alloc() {
    std::unordered_map<DOFManagerType, std::string> types{
        {_dmt_default, "default"}, {_dmt_petsc, "petsc"}};

    return DOFManagerFactory::getInstance().allocate(types[T::value], *mesh,
                                                     "dof_manager", 0);
  }

protected:
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
  dof_manager->registerDOFs("dofs1", *this->dof1, _dst_generic);

  
  
}

TYPED_TEST(DOFManagerFixture, RegisterNodalDOF1) {
  auto dof_manager = this->alloc();
  dof_manager->registerDOFs("dofs1", *this->dof1, _dst_nodal);
}


TYPED_TEST(DOFManagerFixture, RegisterGenericDOF2) {
  auto dof_manager = this->alloc();
  dof_manager->registerDOFs("dofs1", *this->dof1, _dst_generic);
  dof_manager->registerDOFs("dofs2", *this->dof2, _dst_generic);
}

TYPED_TEST(DOFManagerFixture, RegisterNodalDOF2) {
  auto dof_manager = this->alloc();
  dof_manager->registerDOFs("dofs1", *this->dof1, _dst_nodal);
  dof_manager->registerDOFs("dofs2", *this->dof2, _dst_nodal);
}

TYPED_TEST(DOFManagerFixture, RegisterMixedDOF) {
  auto dof_manager = this->alloc();
  dof_manager->registerDOFs("dofs1", *this->dof1, _dst_nodal);
  dof_manager->registerDOFs("dofs2", *this->dof2, _dst_generic);
}
