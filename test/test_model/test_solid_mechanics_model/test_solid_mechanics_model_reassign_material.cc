/**
 * @file   test_solid_mechanics_model_reassign_material.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Mon Feb 10 2014
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  test the function reassign material
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
#include "aka_grid_dynamic.hh"
#include "communicator.hh"
#include "material.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

class StraightInterfaceMaterialSelector : public MaterialSelector {
public:
  StraightInterfaceMaterialSelector(SolidMechanicsModel & model,
                                    UInt horizontal, Real & pos_interface,
                                    const std::string & mat_1_material,
                                    const std::string & mat_2_material)
      : model(model), horizontal(horizontal), pos_interface(pos_interface),
        mat_1_material(mat_1_material), mat_2_material(mat_2_material) {
    Mesh & mesh = model.getMesh();
    UInt spatial_dimension = mesh.getSpatialDimension();

    /// store barycenters of all elements
    barycenters.initialize(mesh, _spatial_dimension = spatial_dimension,
                           _nb_component = spatial_dimension,
                           _with_nb_element = true);

    for_each_element(mesh, [&](auto && el) {
      Vector<Real> bary(barycenters.get(el));
      mesh.getBarycenter(el, bary);
    });
  }

  void setMaterials() {
    mat_ids[0] = model.getMaterialIndex(mat_1_material);
    mat_ids[1] = model.getMaterialIndex(mat_2_material);
  }

  UInt operator()(const Element & elem) override {
    if (not materials_set) {
      setMaterials();
    }
    const Vector<Real> bary = barycenters.get(elem);
    /// check for a given element on which side of the material interface plane
    /// the bary center lies and assign corresponding material
    if (bary(horizontal) < pos_interface) {
      return mat_ids[0];
    }
    return mat_ids[1];
  }

  bool isConditonVerified() {
    /// check if material has been (re)-assigned correctly
    auto & mesh = model.getMesh();
    auto spatial_dimension = mesh.getSpatialDimension();
    for (const auto & type : mesh.elementTypes(spatial_dimension)) {
      auto & mat_indexes = model.getMaterialByElement(type);
      for (auto && data :
           enumerate(make_view(barycenters(type), spatial_dimension))) {
        auto elem = std::get<0>(data);
        auto & bary = std::get<1>(data);
        /// compare element_index_by material to material index that should be
        /// assigned due to the geometry of the interface
        UInt mat_index;
        if (bary(horizontal) < pos_interface) {
          mat_index = mat_ids[0];
        } else {
          mat_index = mat_ids[1];
        }

        if (mat_indexes(elem) != mat_index) {
          /// wrong material index, make test fail
          return false;
        }
      }
    }
    return true;
  }

  void moveInterface(Real & pos_new, UInt horizontal_new) {
    /// update position and orientation of material interface plane
    pos_interface = pos_new;
    horizontal = horizontal_new;
    model.reassignMaterial();
  }

protected:
  SolidMechanicsModel & model;
  ElementTypeMapArray<Real> barycenters;
  std::array<UInt, 2> mat_ids;
  UInt horizontal;
  Real pos_interface;
  bool materials_set{false};
  std::string mat_1_material;
  std::string mat_2_material;
};

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {

  bool test_passed;

  debug::setDebugLevel(dblWarning);
  initialize("two_materials.dat", argc, argv);

  /// specify position and orientation of material interface plane
  Real pos_interface = 0.;

  UInt spatial_dimension = 3;

  const auto & comm = Communicator::getStaticCommunicator();
  Int prank = comm.whoAmI();

  Mesh mesh(spatial_dimension);

  if (prank == 0) {
    mesh.read("cube_two_materials.msh");
  }
  mesh.distribute();

  /// model creation
  SolidMechanicsModel model(mesh);

  /// assign the two different materials using the
  /// StraightInterfaceMaterialSelector

  auto && mat_selector = std::make_shared<StraightInterfaceMaterialSelector>(
      model, _x, pos_interface, "mat_1", "mat_2");

  model.setMaterialSelector(mat_selector);
  model.initFull(_analysis_method = _static);
  MeshUtils::buildFacets(mesh);

  /// check if different materials have been assigned correctly
  test_passed = mat_selector->isConditonVerified();
  if (not test_passed) {
    AKANTU_ERROR("materials not correctly assigned");
    return EXIT_FAILURE;
  }

  model.addDumpField("material_index");
  /// change orientation of material interface plane

  model.dump();
  mat_selector->moveInterface(pos_interface, _y);
  model.dump();

  /// test if material has been reassigned correctly
  test_passed = mat_selector->isConditonVerified();
  if (not test_passed) {
    AKANTU_ERROR("materials not correctly reassigned");
    return EXIT_FAILURE;
  }

  finalize();

  if (prank == 0)
    std::cout << "OK: test passed!" << std::endl;

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
