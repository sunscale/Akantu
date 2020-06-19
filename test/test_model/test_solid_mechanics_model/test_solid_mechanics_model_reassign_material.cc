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
                                    const std::string & mat_1_material,
                                    const std::string & mat_2_material,
                                    bool & horizontal, Real & pos_interface)
      : model(model), mat_1_material(mat_1_material),
        mat_2_material(mat_2_material), horizontal(horizontal),
        pos_interface(pos_interface) {
    Mesh & mesh = model.getMesh();
    UInt spatial_dimension = mesh.getSpatialDimension();

    /// store barycenters of all elements
    barycenters.initialize(mesh, _spatial_dimension = spatial_dimension,
                           _nb_component = spatial_dimension);

    for (auto ghost_type : ghost_types) {
      Element e;
      e.ghost_type = ghost_type;

      for (auto & type : mesh.elementTypes(spatial_dimension, ghost_type)) {
        UInt nb_element = mesh.getNbElement(type, ghost_type);
        e.type = type;
        Array<Real> & barycenter = barycenters(type, ghost_type);
        barycenter.resize(nb_element);

        Array<Real>::iterator<Vector<Real>> bary_it =
            barycenter.begin(spatial_dimension);
        for (UInt elem = 0; elem < nb_element; ++elem) {
          e.element = elem;
          mesh.getBarycenter(e, *bary_it);
          ++bary_it;
        }
      }
    }
  }

  UInt operator()(const Element & elem) {
    UInt spatial_dimension = model.getSpatialDimension();
    const Vector<Real> & bary = barycenters(elem.type, elem.ghost_type)
                                    .begin(spatial_dimension)[elem.element];

    /// check for a given element on which side of the material interface plane
    /// the bary center lies and assign corresponding material
    if (bary(!horizontal) < pos_interface) {
      return model.getMaterialIndex(mat_1_material);
      ;
    }
    return model.getMaterialIndex(mat_2_material);
    ;
  }

  bool isConditonVerified() {

    /// check if material has been (re)-assigned correctly
    Mesh & mesh = model.getMesh();
    UInt spatial_dimension = mesh.getSpatialDimension();
    GhostType ghost_type = _not_ghost;

    for (auto & type : mesh.elementTypes(spatial_dimension, ghost_type)) {
      Array<UInt> & mat_indexes = model.getMaterialByElement(type, ghost_type);
      UInt nb_element = mesh.getNbElement(type, ghost_type);
      Array<Real>::iterator<Vector<Real>> bary =
          barycenters(type, ghost_type).begin(spatial_dimension);
      for (UInt elem = 0; elem < nb_element; ++elem, ++bary) {
        /// compare element_index_by material to material index that should be
        /// assigned due to the geometry of the interface
        UInt mat_index;
        if ((*bary)(!horizontal) < pos_interface)
          mat_index = model.getMaterialIndex(mat_1_material);
        else
          mat_index = model.getMaterialIndex(mat_2_material);

        if (mat_indexes(elem) != mat_index)
          /// wrong material index, make test fail
          return false;
      }
    }
    return true;
  }

  void moveInterface(Real & pos_new, bool & horizontal_new) {
    /// update position and orientation of material interface plane
    pos_interface = pos_new;
    horizontal = horizontal_new;
    model.reassignMaterial();
  }

protected:
  SolidMechanicsModel & model;
  ElementTypeMapArray<Real> barycenters;
  std::string mat_1_material;
  std::string mat_2_material;
  bool horizontal;
  Real pos_interface;
};

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {

  bool test_passed;

  debug::setDebugLevel(dblWarning);
  initialize("two_materials.dat", argc, argv);

  /// specify position and orientation of material interface plane
  bool horizontal = true;
  Real pos_interface = 0.;

  UInt spatial_dimension = 3;

  const auto & comm = Communicator::getStaticCommunicator();
  Int prank = comm.whoAmI();

  Mesh mesh(spatial_dimension);

  if (prank == 0)
    mesh.read("cube_two_materials.msh");
  mesh.distribute();

  /// model creation
  SolidMechanicsModel model(mesh);

  /// assign the two different materials using the
  /// StraightInterfaceMaterialSelector

  auto && mat_selector = std::make_shared<StraightInterfaceMaterialSelector>(
      model, "mat_1", "mat_2", horizontal, pos_interface);

  model.setMaterialSelector(mat_selector);
  model.initFull(_analysis_method = _static);
  MeshUtils::buildFacets(mesh);

  /// check if different materials have been assigned correctly
  test_passed = mat_selector->isConditonVerified();
  if (!test_passed) {
    AKANTU_ERROR("materials not correctly assigned");
    return EXIT_FAILURE;
  }

  /// change orientation of material interface plane
  horizontal = false;
  mat_selector->moveInterface(pos_interface, horizontal);

  // model.dump();

  /// test if material has been reassigned correctly
  test_passed = mat_selector->isConditonVerified();
  if (!test_passed) {
    AKANTU_ERROR("materials not correctly reassigned");
    return EXIT_FAILURE;
  }

  finalize();

  if (prank == 0)
    std::cout << "OK: test passed!" << std::endl;

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
