/**
 * @file   test_.cc
 *
 * @author Clement Roux-Langlois <clement.roux@epfl.ch>
 *
 * @date creation: Fri Mar 13 2015
 * @date last modification: Tue june 16 2015
 *
 * @brief  Tests the interface mesh generation
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "mesh_sphere_intersector.hh"
#include "mesh_igfem_spherical_growing_gel.hh"

#include "dumper_paraview.hh"
#include "dumpable_inline_impl.hh"

/* -------------------------------------------------------------------------- */

using namespace akantu;

typedef Spherical SK;

/* -------------------------------------------------------------------------- */

int main (int argc, char * argv[]) {
  initialize("", argc, argv);
  debug::setDebugLevel(dblError);

  Math::setTolerance(1e-12);

  Mesh mesh(2);
  mesh.read("test_geometry_triangle.msh");
  //mesh.read("mesh.msh");

  // Spherical kernel testing the addition of nodes
  SK::Sphere_3 sphere(SK::Point_3(1, 0, 0), 1.2*1.2); // 0.52); //
  MeshSphereIntersector<2, _triangle_3> intersector_sphere3(mesh);
  MeshSphereIntersector<2, _igfem_triangle_4> intersector_sphere4(mesh);
  MeshSphereIntersector<2, _igfem_triangle_5> intersector_sphere5(mesh);

  std::list<SK::Sphere_3> sphere_list;
  sphere_list.push_back(sphere);

  DumperParaview dumper_igfem("mesh_igfem");
  dumper_igfem.registerMesh(mesh, 2, _not_ghost, _ek_igfem);
  DumperParaview dumper_regular("mesh_regular");
  dumper_regular.registerMesh(mesh, 2, _not_ghost, _ek_regular);
  //dumper_igfem.dump();
  dumper_regular.dump();
  
  intersector_sphere3.buildResultFromQueryList(sphere_list);
  intersector_sphere4.buildResultFromQueryList(sphere_list); 
  intersector_sphere5.buildResultFromQueryList(sphere_list);;


  dumper_regular.dump(); 
  //dumper_igfem.getDumper().setMode(iohelper::TEXT);
  dumper_igfem.dump();

  SK::Sphere_3 sphere2(SK::Point_3(1, 0, 0), 0.4999999999);
  sphere_list.push_back(sphere2);
  intersector_sphere3.constructData();
  intersector_sphere4.constructData();
  intersector_sphere5.constructData();
  intersector_sphere3.removeAdditionnalNodes();
  
  intersector_sphere3.buildResultFromQueryList(sphere_list);
  intersector_sphere4.buildResultFromQueryList(sphere_list);
  intersector_sphere5.buildResultFromQueryList(sphere_list);

  dumper_igfem.dump();
  dumper_regular.dump();
  
  UInt nb_tri3 = mesh.getConnectivity(_triangle_3).getSize();
  UInt nb_tri4 = mesh.getConnectivity(_igfem_triangle_4).getSize();
  UInt nb_tri5 = mesh.getConnectivity(_igfem_triangle_5).getSize();
  std::cout << "final mesh with " << nb_tri3 << " _triangle_3, and " << nb_tri4
	    << " _igfem_triangle_4, and " << nb_tri5 << " _igfem_triangle_5"<< std::endl; 
  if ( (nb_tri3 != 0) || (nb_tri4 != 1) || (nb_tri5 != 1)){ 
    std::cout << "final mesh with " << nb_tri3 << " _triangle_3, and " << nb_tri4
	      << " _igfem_triangle_4, and " << nb_tri5 << " _igfem_triangle_5"<< std::endl; 
    return EXIT_FAILURE;
  }

  /// test of MeshIgfemSphericalGrowingGel
  Mesh mesh_gel(2, "mesh_gel");
  mesh_gel.read("test_geometry_triangle.msh");
  MeshIgfemSphericalGrowingGel<2> gel_intersector(mesh_gel);

  SK::Sphere_3 sphere_gel(SK::Point_3(1, 0, 0), 0.4999999999/4);
  std::list<SK::Sphere_3> sphere_list_gel;
  sphere_list_gel.push_back(sphere_gel);

  DumperParaview dumper_gel_igfem("mesh_gel_igfem");
  dumper_gel_igfem.registerMesh(mesh_gel, 2, _not_ghost, _ek_igfem);
  DumperParaview dumper_gel_regular("mesh_gel_regular");
  dumper_gel_regular.registerMesh(mesh_gel, 2, _not_ghost, _ek_regular);
  dumper_gel_regular.dump();

  gel_intersector.buildResultFromQueryList(sphere_list_gel);
  dumper_gel_igfem.dump();
  dumper_gel_regular.dump();

  UInt nb_tri3_gel = mesh_gel.getConnectivity(_triangle_3).getSize();
  UInt nb_tri4_gel = mesh_gel.getConnectivity(_igfem_triangle_4).getSize();
  UInt nb_tri5_gel = mesh_gel.getConnectivity(_igfem_triangle_5).getSize();
  if ( (nb_tri3_gel != 1) || (nb_tri4_gel != 0) || (nb_tri5_gel != 1)){ 
    std::cout << "interm. mesh_gel with " << nb_tri3_gel << " _triangle_3, and " << nb_tri4_gel
	      << " _igfem_triangle_4, and " << nb_tri5_gel << " _igfem_triangle_5"<< std::endl; 
    return EXIT_FAILURE;
  }

  gel_intersector.buildResultFromQueryList(sphere_list_gel, 1.5);
  dumper_gel_igfem.dump();
  dumper_gel_regular.dump();
  
  gel_intersector.buildResultFromQueryList(sphere_list_gel, 2);
  dumper_gel_igfem.dump();
  dumper_gel_regular.dump();

  gel_intersector.buildResultFromQueryList(sphere_list_gel,2.1);
  dumper_gel_igfem.dump();
  dumper_gel_regular.dump();
  
  gel_intersector.buildResultFromQueryList(sphere_list_gel,3.5);
  dumper_gel_igfem.dump();
  dumper_gel_regular.dump();

  nb_tri3_gel = mesh_gel.getConnectivity(_triangle_3).getSize();
  nb_tri4_gel = mesh_gel.getConnectivity(_igfem_triangle_4).getSize();
  nb_tri5_gel = mesh_gel.getConnectivity(_igfem_triangle_5).getSize();
  if ( (nb_tri3_gel != 1) || (nb_tri4_gel != 0) || (nb_tri5_gel != 1)){ 
    std::cout << "final mesh with " << nb_tri3_gel << " _triangle_3, and " << nb_tri4_gel
	      << " _igfem_triangle_4, and " << nb_tri5_gel << " _igfem_triangle_5"<< std::endl; 
    return EXIT_FAILURE;
    }

  finalize();
  return EXIT_SUCCESS;
}


