/**
 * @file   test_dcb_2d.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Sep 19 14:28:27 2014
 *
 * @brief  2D DCB test to verify the convergence to a same solution
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
#include <limits>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "material.hh"
#include "material_cohesive.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[]) {

  std::cout << " ./test_dcb_2d <final opening> <nb steps> [ <bool_dump>=true  <paraview_folder>=test_dcb_2d ] " << std::endl;

  //debug::setDebugLevel(dblWarning);
  initialize("input_test.dat", argc, argv);
  
  // Math::setTolerance(1.e-13);
  Real tolerance = Math::getTolerance();
  
  const UInt spatial_dimension = 2;
  const Real final_opening =std::atof(argv[1]);
  const UInt step = std::atoi(argv[2]);
  bool bool_dump = true;
  std::string simulation_name = "test_dcb_2d";

  if (argc > 3) bool_dump = std::atoi(argv[3]);
  if (argc > 4) simulation_name = argv[4]; 

  std::cout << "final opening = " << final_opening << " nb_steps = " << step; 
  if (bool_dump) std::cout << " paraview_folder: paraview/" << simulation_name << std::endl;
  else std::cout << std::endl;

  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();
  akantu::MeshPartition * partition = NULL;

  if(prank==0){

    mesh.read("mesh_dcb_2d.msh");

    //CohesiveElementInserter inserter(mesh);
    //inserter.setLimit(_y, -1e-8, 1e-8);
    //inserter.insertIntrinsicElements();
    
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);

  }

  SolidMechanicsModelCohesive model(mesh);

  model.initParallel(partition);

  model.initFull(SolidMechanicsModelCohesiveOptions(_static));

  model.limitInsertion(_y, -1e-8, 1e-8);
  model.insertIntrinsicElements();

  mesh.createGroupsFromMeshData<std::string>("physical_names");
  
  const Array<Real> & nodes = mesh.getNodes();
  Array<UInt> top_boundary_nodes, bot_boundary_nodes;
  Array<bool> & blockedDOFs = model.getBlockedDOFs();
  Array<Real> & displacement = model.getDisplacement();
  UInt nb_nodes = mesh.getNbNodes();

 for (UInt i = 0; i < nb_nodes; ++i) {
     
   if(std::abs(nodes(i,0)-1)<tolerance){
	
     if((nodes(i,1)>0.)&&((nodes(i,1)<0.02))){
       
       top_boundary_nodes.push_back(i);
       blockedDOFs(i,1) = true;
       std::cout << "+"<< std::endl;

     } else if ((nodes(i,1)<0.)&&((nodes(i,1)>-0.02))){
       
       bot_boundary_nodes.push_back(i);
       blockedDOFs(i,1) = true;
       std::cout << "-"<< std::endl;
     }
   }

   if (std::abs(nodes(i,0)) < tolerance) {
  
       blockedDOFs(i,0) = true;
       blockedDOFs(i,1) = true;
   }
 }
 
 model.synchronizeBoundaries();
 model.updateResidual();

 std::cout << mesh.getNbElement(_cohesive_2d_6) << std::endl;

 if (bool_dump) {

   std::stringstream paraview_folder;
   
   paraview_folder << "paraview"
		   <<"/"
		   << simulation_name
		   << "/";
 
   model.setDirectory(paraview_folder.str());
   model.setBaseName("bulk");
   model.addDumpFieldVector("displacement");
   model.addDumpField("stress");
   model.addDumpField("partitions");
   //model.addDumpField("strain");
   model.dump();

   model.setDirectoryToDumper("cohesive elements", "test_unique");
   model.setBaseNameToDumper("cohesive elements", "one_cohesive_element");
   model.addDumpFieldVectorToDumper("cohesive elements", "displacement");
   model.dump("cohesive elements");
 }

 model.assembleStiffnessMatrix();

 Real opening = final_opening/step;

 std::ofstream node_f;
 node_f.open("/home/fabarras/clement_outputs/node_coords.out");
 std::ofstream displ_f;
 displ_f.open("/home/fabarras/clement_outputs/displ.out");

 for (UInt n = 0; n < nb_nodes; ++n) {
   for (UInt d = 0; d < spatial_dimension; ++d) {
       
     node_f << nodes(n,d) << " ";        
   }
   node_f << std::endl;
 }


 for (UInt stp = 0; stp <step; ++stp) {
   std::cout << stp << std::endl;
   for (UInt i = 0; i < top_boundary_nodes.getSize(); ++i) {

     displacement(top_boundary_nodes(i),1) += opening;
   }

   for (UInt i = 0; i < bot_boundary_nodes.getSize(); ++i) {

     displacement(bot_boundary_nodes(i),1) -= opening; 
   }

   Real error;

   model.solveStep<_scm_newton_raphson_tangent, _scc_increment>(1e-8, error, 500);

   std::cout << "Error after convergence: " << error << std::endl;

   if (bool_dump){
     model.dump();
     model.dump("cohesive elements");
   }

   if (stp%5 == 0) {
     for (UInt n = 0; n < nb_nodes; ++n) {
       for (UInt d = 0; d < spatial_dimension; ++d) {
       
	 displ_f << displacement(n,d) << " ";        
       }
       }
     std::cout << stp << std::endl;
   }
 }

 node_f.close();
 displ_f.close();

 finalize();
 
 return EXIT_SUCCESS;
}

 
