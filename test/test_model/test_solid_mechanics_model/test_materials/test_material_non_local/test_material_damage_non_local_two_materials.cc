/**
 * @file   test_material_damage_non_local_two_materials.cc
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Fri Feb 20 15:30:27 2015
 *
 * @brief  test to check if two non-local materials can be used together
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
#include <iostream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  initialize("two_materials_non_local.dat", argc, argv);
  debug::setDebugLevel(akantu::dblWarning);

  const UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  akantu::MeshPartition * partition = NULL;

  StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  if(prank == 0) {

    mesh.read("two_materials.msh");

    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);

    partition->partitionate(psize);
  }

  SolidMechanicsModel model(mesh);
  model.initParallel(partition);
  delete partition;

  model.setBaseName("damage_in_ASR_sample");
  model.addDumpField("partitions");

  model.dump();
  /// assign the materials
  MeshDataMaterialSelector<UInt> * mat_selector;
  mat_selector = new MeshDataMaterialSelector<UInt>("tag_0", model, 1);

  model.setMaterialSelector(*mat_selector);
  model.initFull(SolidMechanicsModelOptions(_static));

  /// displaying the material data
  if (prank == 0) {
    for (UInt m = 0; m  < model.getNbMaterials(); ++m)
      std::cout << model.getMaterial(m) << std::endl;
  }


  //  model.getFEEngine().getMesh().initElementTypeMapArray(quadrature_points_volumes, 1, 0);
  // const MaterialNonLocal<spatial_dimension, BaseWeightFunction> & mat =
  //   dynamic_cast<const MaterialNonLocal<spatial_dimension, BaseWeightFunction> &>(model.getMaterial(0));
  // //  mat.computeQuadraturePointsNeighborhoudVolumes(quadrature_points_volumes);
  // Real radius = mat.getRadius();

  // UInt nb_element  = mesh.getNbElement(TYPE);
  // UInt nb_tot_quad = model.getFEEngine().getNbQuadraturePoints(TYPE) * nb_element;

  // std::cout << mat << std::endl;

  // Array<Real> quads(0, spatial_dimension);
  // quads.resize(nb_tot_quad);

  // model.getFEEngine().interpolateOnQuadraturePoints(mesh.getNodes(),
  // 					       quads, spatial_dimension,
  // 					       TYPE);

  // Array<Real>::iterator< Vector<Real> > first_quad_1 = quads.begin(spatial_dimension);
  // Array<Real>::iterator< Vector<Real> > last_quad_1 = quads.end(spatial_dimension);

  // std::ofstream pout;
  // pout.open("bf_pairs");
  // UInt q1 = 0;

  // Real R = mat.getRadius();

  // for(;first_quad_1 != last_quad_1; ++first_quad_1, ++q1) {
  //   Array<Real>::iterator< Vector<Real> > first_quad_2 = quads.begin(spatial_dimension);
  //   //Array<Real>::iterator< Vector<Real> > last_quad_2 = quads.end(spatial_dimension);
  //   UInt q2 = 0;
  //   for(;first_quad_2 != last_quad_1; ++first_quad_2, ++q2) {
  //     Real d = first_quad_2->distance(*first_quad_1);
  //     if(d <= radius) {
  // 	Real alpha = (1 - d*d/(R*R));
  // 	alpha = alpha*alpha;
  // 	pout << q1 << " " << q2 << " " << alpha << std::endl;
  //     }
  //   }
  // }
  // pout.close();

  // mat.savePairs("cl_pairs");

  // ElementTypeMapArray<Real> constant("constant_value", "test");
  // mesh.initElementTypeMapArray(constant, 1, 0);
  // Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  // Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);
  // for(; it != last_type; ++it) {
  //   UInt nb_quadrature_points = model.getFEEngine().getNbQuadraturePoints(*it);
  //   UInt _nb_element = mesh.getNbElement(*it);

  //   Array<Real> & constant_vect = constant(*it);
  //   constant_vect.resize(_nb_element * nb_quadrature_points);

  //   std::fill_n(constant_vect.storage(), nb_quadrature_points * _nb_element, 1.);
  // }

  // ElementTypeMapArray<Real> constant_avg("constant_value_avg", "test");
  // mesh.initElementTypeMapArray(constant_avg, 1, 0);

  // mat.weightedAvergageOnNeighbours(constant, constant_avg, 1);

  // debug::setDebugLevel(akantu::dblTest);
  // std::cout << constant(TYPE) << std::endl;
  // std::cout << constant_avg(TYPE) << std::endl;
  // debug::setDebugLevel(akantu::dblWarning);

  finalize();

  return EXIT_SUCCESS;
}

