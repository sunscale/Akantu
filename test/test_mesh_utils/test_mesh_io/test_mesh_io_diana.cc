/**
 * @file   test_mesh_io_diana.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Alodie Schneuwly <alodie.schneuwly@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed May 04 18:17:20 2011
 *
 * @brief  test reading mesh diana
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


/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_diana.hh"

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"
#endif //AKANTU_USE_IOHELPER

/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  int dim = 3;
  const ElementType element_type = _tetrahedron_4;
  const UInt paraview_type = iohelper::TETRA1;

  akantu::MeshIODiana mesh_io;
  akantu::Mesh mesh(3);

  mesh_io.read("./dam.dat", mesh);

  std::cout << mesh << std::endl;

  UInt nb_nodes = mesh.getNbNodes();
  UInt nb_elements = mesh.getNbElement(element_type);

  std::vector<std::string> node_groups    = mesh_io.getNodeGroupsNames();
  std::vector<std::string> element_groups = mesh_io.getElementGroupsNames();

#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
  iohelper::DumperParaview dumper;
  dumper.SetMode(iohelper::TEXT);
  dumper.SetPoints(mesh.getNodes().storage(), dim, nb_nodes, "dam_diana");
  dumper.SetConnectivity((int *)mesh.getConnectivity(element_type).storage(), paraview_type, nb_elements, iohelper::C_MODE);

  UInt i = 0;

  Real * nodes_grps[node_groups.size()];
  std::cout << "Nb node groups : " << node_groups.size() << std::endl;
  std::vector<std::string>::iterator it_nodes;
  for(it_nodes = node_groups.begin(); it_nodes != node_groups.end(); ++it_nodes) {
    nodes_grps[i] = new Real[nb_nodes];
    std::fill_n(nodes_grps[i], nb_nodes, 0);
    const Array<UInt> & group = mesh_io.getNodeGroup(*it_nodes);
    for (UInt n = 0; n < group.getSize(); ++n) {
      nodes_grps[i][group(n)] = 1.;
    }

    dumper.AddNodeDataField(nodes_grps[i], 1, (*it_nodes).c_str());

    std::cout << " " << *it_nodes;
    i++;
  }
  std::cout << std::endl;

  i = 0;
  Real * elements_grps[element_groups.size()];
  std::cout << "Nb element groups : " << element_groups.size() << std::endl;
  std::vector<std::string>::iterator it_elements;
  for(it_elements = element_groups.begin(); it_elements != element_groups.end(); ++it_elements) {
    elements_grps[i] = new Real[nb_elements];
    std::fill_n(elements_grps[i], nb_elements, 0);
    const std::vector<Element> & group = mesh_io.getElementGroup(*it_elements);
    for (UInt n = 0; n < group.size(); ++n) {
      elements_grps[i][group[n].element] = 1.;
    }

    dumper.AddElemDataField(elements_grps[i], 1, (*it_elements).c_str());

    std::cout << " " << *it_elements;
    i++;
  }
  std::cout << std::endl;


  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  return EXIT_SUCCESS;
}

