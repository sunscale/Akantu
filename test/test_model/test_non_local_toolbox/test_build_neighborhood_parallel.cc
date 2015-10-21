/**
 * @file   test_build_neighborhood_parallel.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Sun Oct 11 11:20:23 2015
 *
 * @brief  test in parallel for the class NonLocalNeighborhood
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
#include "solid_mechanics_model.hh"
#include "test_material_damage.hh"
#include "non_local_neighborhood_base.hh"
#include "dumper_paraview.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  akantu::initialize("material_parallel_test.dat", argc, argv);

  StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  // some configuration variables
  const UInt spatial_dimension = 2;

  // mesh creation and read
  Mesh mesh(spatial_dimension);
  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {

    mesh.read("fine_mesh.msh");


    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);

    partition->partitionate(psize);
  }

  /// model creation
  SolidMechanicsModel  model(mesh);
  model.initParallel(partition);
  delete partition;


  /// dump the ghost elements before the non-local part is intialized
  DumperParaview dumper_ghost("ghost_elements");
  dumper_ghost.registerMesh(mesh, spatial_dimension, _ghost);
  if(psize > 1) {              
    dumper_ghost.dump();
  }

  /// creation of material selector
  MeshDataMaterialSelector<std::string> * mat_selector;
  mat_selector = new MeshDataMaterialSelector<std::string>("physical_names", model);
  model.setMaterialSelector(*mat_selector);

  /// dump material index in paraview
  model.addDumpField("partitions");
  model.dump();

  /// model initialization changed to use our material
  model.initFull(SolidMechanicsModelOptions(_static, true));
  model.registerNewCustomMaterials< TestMaterialDamage<spatial_dimension> >("test_material");
  model.initMaterials();
  /// dump the ghost elements after ghosts for non-local have been added
  if(psize > 1)
    dumper_ghost.dump();

  model.addDumpField("grad_u");
  model.addDumpField("grad_u non local");
  model.addDumpField("material_index");

 /// apply constant strain field everywhere in the plate
  Matrix<Real> applied_strain(spatial_dimension, spatial_dimension);
  applied_strain.clear();
  for (UInt i = 0; i < spatial_dimension; ++i)
    applied_strain(i,i) = 2.;

  ElementType element_type = _triangle_3;
  GhostType ghost_type = _not_ghost;
  /// apply constant grad_u field in all elements
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    Material & mat = model.getMaterial(m);
    Array<Real> & grad_u = const_cast<Array<Real> &> (mat.getInternal<Real>("grad_u")(element_type, ghost_type));
    Array<Real>::iterator< Matrix<Real> > grad_u_it = grad_u.begin(spatial_dimension, spatial_dimension);
    Array<Real>::iterator< Matrix<Real> > grad_u_end = grad_u.end(spatial_dimension, spatial_dimension);
    /// apply different strain in the first element on each partition
    if (grad_u_it != grad_u_end) {
      (*grad_u_it) += (2. *applied_strain);
      ++grad_u_it;
    }
    for (; grad_u_it != grad_u_end; ++grad_u_it) 
      (*grad_u_it) += applied_strain;
  }

  /// compute the non-local strains
  model.getNonLocalManager().computeAllNonLocalStresses();
  model.dump();


  /// print results to screen for validation
  // std::ifstream quad_pairs;
  // quad_pairs.open("quadrature_pairs.0");
  // std::string current_line;
  // while(getline(quad_pairs, current_line))
  //   std::cout << current_line << std::endl;
  // quad_pairs.close();
  // std::ifstream neighborhoods;
  // neighborhoods.open("neighborhoods.0");
  // while(getline(neighborhoods, current_line))
  //   std::cout << current_line << std::endl;
  // neighborhoods.close();

  finalize();
  
  return EXIT_SUCCESS;
}
