/**
 * @file   blocks.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date   Tue Jan 14 09:38:00 2014
 *
 * @brief  Example of two blocks in contact
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


#include <iostream>
#include <chrono>

#include "contact_impl.hh"

using namespace akantu;


using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  
  // set dimension
  static const UInt dim = 2;
  
  typedef std::chrono::high_resolution_clock clock;
  typedef std::chrono::seconds seconds;
  
  typedef SolidMechanicsModel model_type;
  
  typedef Contact <dim, MasterAssignator, SelectResolution <_static, _augmented_lagrangian> >
  contact_type;

  // capture initial time
  clock::time_point t0 = clock::now();
  
  initialize("steel.dat", argc, argv);
  
  // create and read mesh
  Mesh mesh(dim);
  mesh.read("blocks.msh");
  
  // create model
  model_type model(mesh);
  SolidMechanicsModelOptions opt(_static);
  
  // initialize material
  model.initFull(opt);

  // setup paraview dumper
  model.setBaseName("contact");
  model.addDumpFieldVector("displacement");
  
  // create data structure that holds contact data
  contact_type cd(argc, argv, model);
  
  // get areas for slave nodes
  mesh.createGroupsFromMeshData<std::string>("physical_names");
  model.applyBC(BC::Neumann::FromHigherDim(Matrix<Real>::eye(2,1.)), "interface_top");
  Array<Real>& areas = model.getForce();
  
  ElementGroup &eg = mesh.getElementGroup("interface_top");
  for (auto nit = eg.node_begin(); nit != eg.node_end(); ++nit) {
    
    // add slave node
    cd.addSlave(*nit);

    // compute and add area
    Real a = 0.;
    for (UInt i=0; i<dim; ++i)
      a += pow(areas(*nit, i),2.);
    cd.addArea(*nit, sqrt(a));
  }
  
  // set force value to zero
  areas.clear();
  
  // add master surface to find pairs
	cd.searchSurface("interface_bottom");

  // apply boundary conditions
  using BC::Dirichlet::FixedValue;
  model.applyBC(FixedValue(0., _x), "bottom");
  model.applyBC(FixedValue(0., _y), "bottom");
  model.applyBC(FixedValue(0., _x), "top");
  
  Real U = 0.5;
  Real Du = 0.005;
  // loop over command line parameters
  for (int i=0; i<argc; ++i) {
    if (strcmp(argv[i], "-steps") == 0) {
      Real steps = atof(argv[++i]);
      cout<<"-steps = "<<steps<<endl;
      Du = U / steps;
    }
  }
  
  // loop over load increments
  for (Real u = Du; u<=U; u += Du) {
    
    model.applyBC(FixedValue(-u, _y), "top");
    solveContactStep<_generalized_newton>(cd);
  }
  
  clock::time_point t1 = clock::now();
  
  seconds total_s = std::chrono::duration_cast<seconds>(t1 - t0);
  
  cout<<"- Simulation took "<<total_s.count()<<" s"<<endl;
  
  // finalize simulation
  finalize();
  return EXIT_SUCCESS;
}

