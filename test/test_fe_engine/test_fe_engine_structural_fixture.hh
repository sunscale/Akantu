/**
 * @file   test_fe_engine_structural_fixture.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Aug 20 2010
 * @date last modification: Fri Jan 26 2018
 *
 * @brief  test of the fem class
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh_io_msh_struct.hh"
#include "test_fe_engine_fixture.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TEST_FE_ENGINE_STRUCTURAL_FIXTURE_HH__
#define __AKANTU_TEST_FE_ENGINE_STRUCTURAL_FIXTURE_HH__

using namespace akantu;

/// Base class for structural FEEngine tests with structural elements
template <typename type_>
class TestFEMStructuralFixture
    : public TestFEMBaseFixture<type_, ShapeStructural, _ek_structural> {
  using parent = TestFEMBaseFixture<type_, ShapeStructural, _ek_structural>;

public:
  static const UInt ndof = ElementClass<parent::type>::getNbDegreeOfFreedom();

  /// Need to tell the mesh to load structural elements
  void readMesh(std::string file_name) override {
    this->mesh->read(file_name, _miot_gmsh_struct);
  }
};

template <typename type_> const UInt TestFEMStructuralFixture<type_>::ndof;

// using types = gtest_list_t<TestElementTypes>;

// TYPED_TEST_SUITE(TestFEMFixture, types);

#endif /* __AKANTU_TEST_FE_ENGINE_STRUCTURAL_FIXTURE_HH__ */
