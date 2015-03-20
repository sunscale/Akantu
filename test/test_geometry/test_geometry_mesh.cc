/**
 * @file   test_geometry_mesh.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Mar 13 2015
 * @date last modification: Fri Mar 13 2015
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

#include "mesh_geom_container.hh"
#include "geom_helper_functions.hh"

#include <CGAL/Cartesian.h>

#include <iostream>

/* -------------------------------------------------------------------------- */

using namespace akantu;

typedef CGAL::Cartesian<Real> K;

/* -------------------------------------------------------------------------- */

int main (int argc, char * argv[]) {
  debug::setDebugLevel(dblWarning);
  initialize("", argc, argv);

  Math::setTolerance(1e-10);

  Mesh mesh(2);
  mesh.read("mesh.msh");

  MeshGeomContainer container(mesh);
  container.constructData();

  K::Point_3 a(0, 0.25, 0),
             b(1, 0.25, 0),
             c(0.25, 0, 0),
             d(0.25, 1, 0);

  K::Segment_3 h_interface(a, b),
               v_interface(c, d);

  std::list<K::Segment_3> interface_list;
  interface_list.push_back(h_interface);
  interface_list.push_back(v_interface);

  Mesh & interface_mesh = container.meshOfLinearInterfaces(interface_list);

  if (interface_mesh.getNbElement(_segment_2) != 4)
    return EXIT_FAILURE;

  Vector<Real> bary(2);
  Element test;
  test.element = 1;
  test.type = _segment_2;
  
  interface_mesh.getBarycenter(test, bary);

  if (!Math::are_float_equal(bary(0), 0.5) || !Math::are_float_equal(bary(1), 0.25))
    return EXIT_FAILURE;

  finalize();
  return EXIT_SUCCESS;
}


