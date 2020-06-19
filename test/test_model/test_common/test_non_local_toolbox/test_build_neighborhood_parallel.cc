/**
 * @file   test_build_neighborhood_parallel.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sat Sep 26 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  test in parallel for the class NonLocalNeighborhood
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dumper_iohelper_paraview.hh"
#include "non_local_neighborhood_base.hh"
#include "solid_mechanics_model.hh"
#include "test_material.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  akantu::initialize("material_parallel_test.dat", argc, argv);

  const auto & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  // some configuration variables
  const UInt spatial_dimension = 2;

  // mesh creation and read
  Mesh mesh(spatial_dimension);
  if (prank == 0) {
    mesh.read("parallel_test.msh");
  }

  mesh.distribute();

  /// model creation
  SolidMechanicsModel model(mesh);

  /// dump the ghost elements before the non-local part is intialized
  DumperParaview dumper_ghost("ghost_elements");
  dumper_ghost.registerMesh(mesh, spatial_dimension, _ghost);
  if (psize > 1) {
    dumper_ghost.dump();
  }

  /// creation of material selector
  auto && mat_selector =
      std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names",
                                                              model);
  model.setMaterialSelector(mat_selector);

  /// dump material index in paraview
  model.addDumpField("partitions");
  model.dump();

  /// model initialization changed to use our material
  model.initFull();

  /// dump the ghost elements after ghosts for non-local have been added
  if (psize > 1)
    dumper_ghost.dump();

  model.addDumpField("grad_u");
  model.addDumpField("grad_u non local");
  model.addDumpField("material_index");

  /// apply constant strain field everywhere in the plate
  Matrix<Real> applied_strain(spatial_dimension, spatial_dimension);
  applied_strain.clear();
  for (UInt i = 0; i < spatial_dimension; ++i)
    applied_strain(i, i) = 2.;

  ElementType element_type = _triangle_3;
  GhostType ghost_type = _not_ghost;
  /// apply constant grad_u field in all elements
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    auto & mat = model.getMaterial(m);
    auto & grad_u = const_cast<Array<Real> &>(
        mat.getInternal<Real>("grad_u")(element_type, ghost_type));
    auto grad_u_it = grad_u.begin(spatial_dimension, spatial_dimension);
    auto grad_u_end = grad_u.end(spatial_dimension, spatial_dimension);
    for (; grad_u_it != grad_u_end; ++grad_u_it)
      (*grad_u_it) = -1. * applied_strain;
  }

  /// double the strain in the center: find the closed gauss point to the center
  /// compute the quadrature points
  ElementTypeMapReal quad_coords("quad_coords");
  quad_coords.initialize(mesh, _nb_component = spatial_dimension,
                         _spatial_dimension = spatial_dimension,
                         _with_nb_element = true);
  model.getFEEngine().computeIntegrationPointsCoordinates(quad_coords);

  Vector<Real> center(spatial_dimension, 0.);
  Real min_distance = 2;
  IntegrationPoint q_min;
  for (auto type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_regular)) {
    UInt nb_elements = mesh.getNbElement(type, _not_ghost);
    UInt nb_quads = model.getFEEngine().getNbIntegrationPoints(type);
    Array<Real> & coords = quad_coords(type, _not_ghost);
    auto coord_it = coords.begin(spatial_dimension);
    for (UInt e = 0; e < nb_elements; ++e) {
      for (UInt q = 0; q < nb_quads; ++q, ++coord_it) {
        Real dist = center.distance(*coord_it);
        if (dist < min_distance) {
          min_distance = dist;
          q_min.element = e;
          q_min.num_point = q;
          q_min.global_num = nb_elements * nb_quads + q;
          q_min.type = type;
        }
      }
    }
  }

  Real global_min = min_distance;
  comm.allReduce(global_min, SynchronizerOperation::_min);

  if (Math::are_float_equal(global_min, min_distance)) {
    UInt mat_index = model.getMaterialByElement(q_min.type, _not_ghost)
                         .begin()[q_min.element];
    Material & mat = model.getMaterial(mat_index);
    UInt nb_quads = model.getFEEngine().getNbIntegrationPoints(q_min.type);
    UInt local_el_index =
        model.getMaterialLocalNumbering(q_min.type, _not_ghost)
            .begin()[q_min.element];
    UInt local_num = (local_el_index * nb_quads) + q_min.num_point;
    Array<Real> & grad_u = const_cast<Array<Real> &>(
        mat.getInternal<Real>("grad_u")(q_min.type, _not_ghost));
    Array<Real>::iterator<Matrix<Real>> grad_u_it =
        grad_u.begin(spatial_dimension, spatial_dimension);
    grad_u_it += local_num;
    Matrix<Real> & g_u = *grad_u_it;
    g_u += applied_strain;
  }

  /// compute the non-local strains
  model.assembleInternalForces();
  model.dump();

  /// damage the element with higher grad_u completely, so that it is
  /// not taken into account for the averaging
  if (Math::are_float_equal(global_min, min_distance)) {
    UInt mat_index = model.getMaterialByElement(q_min.type, _not_ghost)
                         .begin()[q_min.element];
    Material & mat = model.getMaterial(mat_index);
    UInt nb_quads = model.getFEEngine().getNbIntegrationPoints(q_min.type);
    UInt local_el_index =
        model.getMaterialLocalNumbering(q_min.type, _not_ghost)
            .begin()[q_min.element];
    UInt local_num = (local_el_index * nb_quads) + q_min.num_point;
    Array<Real> & damage = const_cast<Array<Real> &>(
        mat.getInternal<Real>("damage")(q_min.type, _not_ghost));
    Real * dam_ptr = damage.storage();
    dam_ptr += local_num;
    *dam_ptr = 0.9;
  }

  /// compute the non-local strains
  model.assembleInternalForces();
  model.dump();

  finalize();

  return EXIT_SUCCESS;
}
