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

#include "aka_common.hh"

#include "mesh_sphere_intersector.hh"
#include "geom_helper_functions.hh"
#include "mesh_geom_common.hh"

#include "dumper_paraview.hh"

#include <iostream>

/* -------------------------------------------------------------------------- */

using namespace akantu;

typedef Spherical SK;

/* -------------------------------------------------------------------------- */

int main (int argc, char * argv[]) {
  initialize("", argc, argv);
  debug::setDebugLevel(dblError);

  Math::setTolerance(1e-10);

  Mesh mesh(2);
  mesh.read("test_geometry_triangle.msh");

  // Spherical kernel testing the addition of nodes
  SK::Sphere_3 sphere(SK::Point_3(0, 1, 0), 0.2*0.2);
  SK::Sphere_3 sphere2(SK::Point_3(1, 0, 0), 0.4999999999);
  MeshSphereIntersector<2, _triangle_3> intersector_sphere(mesh);
  intersector_sphere.constructData();

  std::list<SK::Sphere_3> sphere_list;
  sphere_list.push_back(sphere);
  sphere_list.push_back(sphere2);

  DumperParaview dumper_igfem("mesh_igfem");
  dumper_igfem.registerMesh(mesh, 2, _not_ghost, _ek_igfem);
  DumperParaview dumper_regular("mesh_regular");
  dumper_regular.registerMesh(mesh, 2, _not_ghost, _ek_regular);
  //dumper_igfem.dump();
  dumper_regular.dump();

  intersector_sphere.buildIgfemMesh(sphere_list);
  dumper_igfem.dump();
  dumper_regular.dump();

  UInt nb_tri3 = mesh.getConnectivity(_triangle_3).getSize();
  UInt nb_tri4 = mesh.getConnectivity(_igfem_triangle_4).getSize();
  UInt nb_tri5 = mesh.getConnectivity(_igfem_triangle_5).getSize();   
  if ( (nb_tri3 != 0) || (nb_tri4 != 1) || (nb_tri5 != 1)){ 
    std::cout << "final mesh with " << nb_tri3 << " _triangle_3, and " << nb_tri4
	      << " _igfem_triangle_4, and " << nb_tri5 << " _igfem_triangle_5"<< std::endl; 
    return EXIT_FAILURE;
  }

  const Array<UInt> new_node_triangle_3 = intersector_sphere.getNewNodePerElem();
  if ( (new_node_triangle_3(0,0) != 1) || (new_node_triangle_3(1,0) != 2)){
    for(UInt k=0; k != new_node_triangle_3.getSize(); ++k){
      std::cout << new_node_triangle_3(k,0) << " new nodes in element " << k << ", node(s): "
		<< new_node_triangle_3(k,1) << ", " << new_node_triangle_3(k,3)
		<< ", on segment(s):" << new_node_triangle_3(k,2) << ", "
		<< new_node_triangle_3(k,4) << std::endl;
    }
    return EXIT_FAILURE;
  }

  finalize();
  return EXIT_SUCCESS;
}


