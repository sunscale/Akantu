/**
 * @file   test_structural_mechanics_model_bernoulli_beam_2_exemple_1_1_y.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Wed Jan 16 2013
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  Test for _bernouilli_beam_2D
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <limits>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "structural_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#define TYPE _bernoulli_beam_2

using namespace akantu;

//Linear load function
static void lin_load(double * position, double * load,
		     __attribute__ ((unused)) Real * normal, __attribute__ ((unused)) UInt surface_id){
  memset(load,0,sizeof(Real)*3);
  if (position[1]<=10){
  load[0]= 6000;
 }
}

int main(int argc, char *argv[]){
  initialize(argc, argv);
  Mesh beams(2);
  debug::setDebugLevel(dblWarning);

  /* -------------------------------------------------------------------------- */
  // Defining the mesh

  UInt nb_nodes=3;
  UInt nb_nodes_1=1;
  UInt nb_nodes_2=nb_nodes-nb_nodes_1 - 1;
  UInt nb_element=nb_nodes-1;

  Array<Real> & nodes = const_cast<Array<Real> &>(beams.getNodes());
  nodes.resize(nb_nodes);

  beams.addConnectivityType(_bernoulli_beam_2);
  Array<UInt> & connectivity = const_cast<Array<UInt> &>(beams.getConnectivity(_bernoulli_beam_2));
  connectivity.resize(nb_element);

  for(UInt i=0; i<nb_nodes; ++i) {
     nodes(i,0)=0;
  }
   for (UInt i = 0; i < nb_nodes_1; ++i) {
     nodes(i,1)=10.*i/((Real)nb_nodes_1);

   }
   nodes(nb_nodes_1,1)=10;

   for (UInt i = 0; i < nb_nodes_2; ++i) {
     nodes(nb_nodes_1 + i + 1,1)=10+ 8.*(i+1)/((Real)nb_nodes_2);

   }

  for(UInt i=0; i<nb_element; ++i) {

    connectivity(i,0)=i;
    connectivity(i,1)=i+1;
  }
  akantu::MeshIOMSH mesh_io;
  mesh_io.write("b_beam_2.msh", beams);

  /* -------------------------------------------------------------------------- */
  // Defining the materials
  akantu::StructuralMechanicsModel model(beams);

  StructuralMaterial mat1;
  mat1.E=3e10;
  mat1.I=0.0025;
  mat1.A=0.01;

  model.addMaterial(mat1);

  StructuralMaterial mat2 ;
  mat2.E=3e10;
  mat2.I=0.00128;
  mat2.A=0.01;

  model.addMaterial(mat2);

  /* -------------------------------------------------------------------------- */
  // Defining the forces
  model.initFull();

  const Real M = -3600; // Momentum at 3

  Array<Real> & forces = model.getForce();
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> & boundary = model.getBlockedDOFs();
  const Array<Real> & N_M  = model.getStress(_bernoulli_beam_2);

  Array<UInt> & element_material = model.getElementMaterial(_bernoulli_beam_2);

  for (UInt i = 0; i < nb_nodes_2; ++i) {
    element_material(i+nb_nodes_1)=1;
  }


  forces(nb_nodes-1,2) += M;
  model.computeForcesFromFunction<_bernoulli_beam_2>(lin_load, akantu::_bft_traction);

  /* -------------------------------------------------------------------------- */
  // Defining the boundary conditions
  boundary(0,0) = true;
  boundary(0,1) = true;
  boundary(0,2) = true;
  boundary(nb_nodes_1,0) = true;
  boundary(nb_nodes-1,0) = true;

  /* -------------------------------------------------------------------------- */
  // Solve
  Real error;

  model.assembleStiffnessMatrix();
  model.getStiffnessMatrix().saveMatrix("Kb.mtx");
  UInt count = 0;

  model.addDumpField("displacememt");
  model.addDumpField("rotation");
  model.addDumpField("force");
  model.addDumpField("momentum");

  do {
    if(count != 0) std::cerr << count << " - " << error << std::endl;
    model.updateResidual();
    model.solve();
    count++;
  } while (!model.testConvergenceIncrement(1e-10, error) && count < 10);
  std::cerr << count << " - " << error << std::endl;

  /* -------------------------------------------------------------------------- */
  // Post-Processing

  model.computeStresses();

  model.getStiffnessMatrix().saveMatrix("Ka.mtx");
  std::cout<< " d1 = " << displacement(nb_nodes_1,2) << std::endl;
  std::cout<< " d2 = " << displacement(nb_nodes-1,2) << std::endl;
  std::cout<< " M1 = " << N_M(0,1) << std::endl;
  std::cout<< " M2 = " << N_M(2*(nb_nodes-2),1) << std::endl;

  model.dump();
}
