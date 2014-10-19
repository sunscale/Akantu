/**
 * @file   tetrahedra.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date   Tue Jan 14 09:38:00 2014
 *
 * @brief  File used to obtain contact results for a simple tetrahedra test
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

#include "contact_impl.hh"

//#include "implicit_contact_manager.hh"

using namespace akantu;

int main(int argc, char *argv[]) {
  
  // set dimension
  static const UInt dim = 3;
  
  typedef SolidMechanicsModel model_type;
//  typedef ContactData<dim,model_type> contact_type;
  typedef Contact <dim, MasterAssignator, SelectResolution <_static, _augmented_lagrangian> >
  contact_type;

  
  // initialize akantu
  initialize("steel.dat", argc, argv);
  
  // create and read mesh
  Mesh mesh(dim);
  mesh.read("tetrahedra.msh");
  
  // create model
  model_type model(mesh);
  SolidMechanicsModelOptions opt(_static);
  
  // initialize model
  model.initFull(opt);
  
  // paraview output
  model.setBaseName("contact");
  model.addDumpFieldVector("displacement");
  
  // create conctact object
  contact_type cd(argc, argv, model);
  
  // add slave node of the tip of the top tetrahedron
  cd.addSlave(4);
  
  // parameters
  cd[Verbose] = true;

  // add search surface
  cd.searchSurface("master");

  // apply boundary conditions
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  using BC::Dirichlet::FixedValue;
//  model.applyBC(FixedValue(0., _y), "top_surface");
  model.applyBC(FixedValue(0., _x), "top_surface");
  model.applyBC(FixedValue(0., _z), "top_surface");
  model.applyBC(FixedValue(0., _x), "bottom_surface");
  model.applyBC(FixedValue(0., _y), "bottom_surface");
  model.applyBC(FixedValue(0., _z), "bottom_surface");

//  model.applyBC(FixedValue(0., _x), "master");
//  model.applyBC(FixedValue(0., _y), "master");
//  model.applyBC(FixedValue(0., _z), "master");

  Real U = 2;
  Real Du = 0.01;
  
  for (Real u = Du; u<=U; u += Du) {
    
//    model.applyBC(FixedValue(u, _x), "top_surface");
    model.applyBC(FixedValue(-u, _y), "top_surface");

    // solve contact step (no need to call solve on the model object)
//    solveContactStep<_uzawa>(cd);
    solveContactStep<_generalized_newton>(cd);
  }
  
  
  // finalize simulation
  finalize();
  return EXIT_SUCCESS;
}
