/**
 * @file   offset_2slaves.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Mon Sep 15 2014
 * @date last modification: Fri Sep 19 2014
 *
 * @brief  This file tests the basic contact implementation using the
 * generalized Newton method
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

#include "contact_impl.hh"

using namespace akantu;

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;

int main(int argc, char *argv[]) {
  
  // set dimension
  static const UInt dim = 2;
  
  // type definitions
  typedef SolidMechanicsModel model_type;
  typedef Contact<dim, MasterAssignator,
  SelectResolution<_static, _augmented_lagrangian> >
  contact_type;
  
  initialize("steel.dat", argc, argv);
  
  // create meshes
  Mesh mesh(dim);
  
  // read meshes
  mesh.read("offset_2slaves.msh");
  
  // create models
  model_type model(mesh);
  
  SolidMechanicsModelOptions opt(_static);
  
  // initialize material
  // initialize material
  model.initFull(opt);
  model.updateCurrentPosition();
  // create data structure that holds contact data
  contact_type cd(argc, argv, model);
  
  // set Paraview output resluts
  model.setBaseName("contact");
  model.addDumpFieldVector("displacement");
  
  mesh.createGroupsFromMeshData<std::string>("physical_names");
  
  cd.addSlave(4);
  cd.addArea(4, 1.0);
  
  // add master surface to find pairs
  cd.searchSurface("Contact");
  
  model.applyBC(BC::Dirichlet::FixedValue(0., _x), "Top");
  model.applyBC(BC::Dirichlet::FixedValue(0., _x), "Bottom");
  model.applyBC(BC::Dirichlet::FixedValue(0., _y), "Bottom");
  
  Real U = 0.5;
  Real Du = 0.01;
  for (Real u = Du; u <= U; u += Du) {
    
    model.applyBC(BC::Dirichlet::FixedValue(-u, _y), "Top");
    
    // solve contact step (no need to call solve on the model object)
    solveContactStep<_generalized_newton>(cd);
  }
  
  cout<<"Force: "<< cd.getForce()<<endl;
  
  // finalize simulation
  finalize();
  return EXIT_SUCCESS;
}
