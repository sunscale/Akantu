/**
 * @file   test_structural_mechanics_model_bernoulli_beam_2_complicated.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Thu Jun 12 2014
 *
 * @brief  A very complicated structure
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
#include "mesh_io_msh_struct.hh"
#include "structural_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#define TYPE _bernoulli_beam_2

using namespace akantu;

//Linear load function
static void lin_load(double * position, double * load,
		     __attribute__ ((unused)) Real * normal, __attribute__ ((unused)) UInt surface_id){
  memset(load,0,sizeof(Real)*3);
  if(position[1]>=0.-Math::getTolerance()) {
    if ((position[0]<=10.)){
      load[1]= -100;
    } else if (position[0]<=20.){
      load[1]= -70;
    }
  }
}

int main(int argc, char *argv[]){

  initialize(argc, argv);
  Mesh beams(2);
  debug::setDebugLevel(dblWarning);

  /* -------------------------------------------------------------------------- */
  // Defining the mesh

  akantu::MeshIOMSHStruct mesh_io;
  mesh_io.read("complicated.msh", beams);

  /* -------------------------------------------------------------------------- */
  // Defining the material

  const akantu::ElementType type = akantu::_bernoulli_beam_2;

  akantu::StructuralMechanicsModel  model(beams);

  StructuralMaterial mat1;
  mat1.E=3e10;
  mat1.I=0.0025;
  mat1.A=0.01;


  model.addMaterial(mat1);

  StructuralMaterial mat2 ;
  mat2.E=3e10;
  mat2.I=0.003125;
  mat2.A=0.01;

  model.addMaterial(mat2);

  /* -------------------------------------------------------------------------- */
  // Defining the forces
  model.initFull();

  UInt nb_element = beams.getNbElement(type);
  for (unsigned int i = 0; i < nb_element; ++i) {
    model.getElementMaterial(type)(i,0) = beams.getData<UInt>("tag_0", type)(i,0) - 1;
  }


  Array<Real> & forces = model.getForce();
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> & boundary = model.getBlockedDOFs();

  forces.clear();
  displacement.clear();

  model.computeForcesFromFunction<_bernoulli_beam_2>(lin_load, akantu::_bft_traction);

  /* -------------------------------------------------------------------------- */
  // Defining the boundary conditions

  boundary(0,0) = true;
  boundary(0,1) = true;
  boundary(3,0) = true;
  boundary(3,1) = true;
  boundary(4,0) = true;
  boundary(4,1) = true;
  boundary(4,2) = true;
  boundary(5,0) = true;
  boundary(5,1) = true;
  boundary(5,2) = true;
  boundary(2,1) = true;
  boundary(2,0) = true;
  boundary(1,1) = true;
  boundary(1,0) = true;


  /* -------------------------------------------------------------------------- */
  // Solve
  Real error;

  model.assembleStiffnessMatrix();
  model.getStiffnessMatrix().saveMatrix("Kb.mtx");
  UInt count = 0;

  model.addDumpFieldVector("displacement");
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
  std::cout<< " x1 = " << displacement(1,2) << std::endl;
  std::cout<< " x2 = " << displacement(2,2) << std::endl;

  model.dump();
}
