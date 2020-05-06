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

#include "mesh_igfem_spherical_growing_gel.hh"
#include "mesh_sphere_intersector.hh"

#include "dumpable_inline_impl.hh"
#include "dumper_paraview.hh"

/* -------------------------------------------------------------------------- */

using namespace akantu;

typedef Spherical SK;

/* -------------------------------------------------------------------------- */

int main(int argc, char * argv[]) {
  initialize("", argc, argv);
  debug::setDebugLevel(dblError);

  Math::setTolerance(1e-12);

  Mesh mesh_gel2(2, "mesh_gel2");
  mesh_gel2.read("test_geometry_triangle.msh");
  MeshIgfemSphericalGrowingGel<2> gel_intersector2(mesh_gel2);
  gel_intersector2.init();

  SK::Sphere_3 sphere_gel2(SK::Point_3(1, 0, 0), 0.4999999999 / 4);
  std::list<SK::Sphere_3> sphere_list_gel2;
  sphere_list_gel2.push_back(sphere_gel2);

  DumperParaview dumper_gel2_igfem("mesh_gel2_igfem");
  dumper_gel2_igfem.registerMesh(mesh_gel2, 2, _not_ghost, _ek_igfem);
  DumperParaview dumper_gel2_regular("mesh_gel2_regular");
  dumper_gel2_regular.registerMesh(mesh_gel2, 2, _not_ghost, _ek_regular);
  dumper_gel2_regular.dump();

  gel_intersector2.buildIGFEMMeshFromSpheres(sphere_list_gel2);
  // gel_intersector2.computeMeshQueryListIntersectionPoint(sphere_list_gel2);
  // gel_intersector2.buildIgfemMesh();
  dumper_gel2_igfem.dump();
  dumper_gel2_regular.dump();
  UInt nb_tri3_gel2 = mesh_gel2.getConnectivity(_triangle_3).getSize();
  UInt nb_tri4_gel2 = mesh_gel2.getConnectivity(_igfem_triangle_4).getSize();
  UInt nb_tri5_gel2 = mesh_gel2.getConnectivity(_igfem_triangle_5).getSize();
  if ((nb_tri3_gel2 != 1) || (nb_tri4_gel2 != 0) || (nb_tri5_gel2 != 1)) {
    std::cout << "interm. mesh_gel with " << nb_tri3_gel2
              << " _triangle_3, and " << nb_tri4_gel2
              << " _igfem_triangle_4, and " << nb_tri5_gel2
              << " _igfem_triangle_5" << std::endl;
    return EXIT_FAILURE;
  }

  gel_intersector2.buildIGFEMMeshFromSpheres(sphere_list_gel2, 1.5);
  // gel_intersector2.computeMeshQueryListIntersectionPoint(sphere_list_gel2,
  // 1.5);
  // gel_intersector2.buildIgfemMesh();
  dumper_gel2_igfem.dump();
  dumper_gel2_regular.dump();
  nb_tri3_gel2 = mesh_gel2.getConnectivity(_triangle_3).getSize();
  nb_tri4_gel2 = mesh_gel2.getConnectivity(_igfem_triangle_4).getSize();
  nb_tri5_gel2 = mesh_gel2.getConnectivity(_igfem_triangle_5).getSize();
  if ((nb_tri3_gel2 != 1) || (nb_tri4_gel2 != 0) || (nb_tri5_gel2 != 1)) {
    std::cout << "final mesh with " << nb_tri3_gel2 << " _triangle_3, and "
              << nb_tri4_gel2 << " _igfem_triangle_4, and " << nb_tri5_gel2
              << " _igfem_triangle_5" << std::endl;
    return EXIT_FAILURE;
  }

  gel_intersector2.buildIGFEMMeshFromSpheres(sphere_list_gel2, 2);
  // gel_intersector2.computeMeshQueryListIntersectionPoint(sphere_list_gel2,
  // 2);
  // gel_intersector2.buildIgfemMesh();
  dumper_gel2_igfem.dump();
  dumper_gel2_regular.dump();
  nb_tri3_gel2 = mesh_gel2.getConnectivity(_triangle_3).getSize();
  nb_tri4_gel2 = mesh_gel2.getConnectivity(_igfem_triangle_4).getSize();
  nb_tri5_gel2 = mesh_gel2.getConnectivity(_igfem_triangle_5).getSize();
  if ((nb_tri3_gel2 != 1) || (nb_tri4_gel2 != 1) || (nb_tri5_gel2 != 0)) {
    std::cout << "final mesh with " << nb_tri3_gel2 << " _triangle_3, and "
              << nb_tri4_gel2 << " _igfem_triangle_4, and " << nb_tri5_gel2
              << " _igfem_triangle_5" << std::endl;
    return EXIT_FAILURE;
  }

  gel_intersector2.buildIGFEMMeshFromSpheres(sphere_list_gel2, 2.1);
  // gel_intersector2.computeMeshQueryListIntersectionPoint(sphere_list_gel2,2.1);
  // gel_intersector2.buildIgfemMesh();
  dumper_gel2_igfem.dump();
  dumper_gel2_regular.dump();
  nb_tri3_gel2 = mesh_gel2.getConnectivity(_triangle_3).getSize();
  nb_tri4_gel2 = mesh_gel2.getConnectivity(_igfem_triangle_4).getSize();
  nb_tri5_gel2 = mesh_gel2.getConnectivity(_igfem_triangle_5).getSize();
  if ((nb_tri3_gel2 != 0) || (nb_tri4_gel2 != 0) || (nb_tri5_gel2 != 2)) {
    std::cout << "final mesh with " << nb_tri3_gel2 << " _triangle_3, and "
              << nb_tri4_gel2 << " _igfem_triangle_4, and " << nb_tri5_gel2
              << " _igfem_triangle_5" << std::endl;
    return EXIT_FAILURE;
  }

  gel_intersector2.buildIGFEMMeshFromSpheres(sphere_list_gel2, 3.5);
  // gel_intersector2.computeMeshQueryListIntersectionPoint(sphere_list_gel,3.5);
  // gel_intersector2.buildIgfemMesh();
  dumper_gel2_igfem.dump();
  dumper_gel2_regular.dump();
  nb_tri3_gel2 = mesh_gel2.getConnectivity(_triangle_3).getSize();
  nb_tri4_gel2 = mesh_gel2.getConnectivity(_igfem_triangle_4).getSize();
  nb_tri5_gel2 = mesh_gel2.getConnectivity(_igfem_triangle_5).getSize();
  if ((nb_tri3_gel2 != 1) || (nb_tri4_gel2 != 0) || (nb_tri5_gel2 != 1)) {
    std::cout << "final mesh with " << nb_tri3_gel2 << " _triangle_3, and "
              << nb_tri4_gel2 << " _igfem_triangle_4, and " << nb_tri5_gel2
              << " _igfem_triangle_5" << std::endl;
    return EXIT_FAILURE;
  }

  finalize();
  return EXIT_SUCCESS;
}
