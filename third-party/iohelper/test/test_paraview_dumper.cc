/**
 * @file   test_paraview_dumper.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Oct 11 2012
 * @date last modification: Tue May 14 2013
 *
 * @brief  paraview dump test
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * IOHelper is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "../src/dumper_paraview.hh"
#include "mesh_io_msh.hh"
/* -------------------------------------------------------------------------- */

int main(int argc, char ** argv){

  if (argc != 3){
    std::cerr << "Usage: program meshfile dimension" << std::endl;
    return EXIT_FAILURE;
  }

  std::string meshfilename(argv[1]);
  int dim = atoi(argv[2]);
 
  std::vector<double> nodes;
  std::vector<int> connectivities;
  int nb_elements;

  MeshIOMSH meshio;
  meshio.read(meshfilename,dim,MeshIOMSH::_msh_tetrahedron_4,nodes,connectivities);
  nb_elements = connectivities.size()/meshio._msh_nodes_per_elem[MeshIOMSH::_msh_tetrahedron_4];

  iohelper::DumperParaview dumper;
  //  dumper.setMode(iohelper::TEXT);
  dumper.setPoints(&nodes[0], dim, nodes.size()/dim, "coordinates2");
  dumper.setConnectivity(&connectivities[0],iohelper::TETRA1, nb_elements, iohelper::C_MODE);
  dumper.addNodeDataField("test",&nodes[0],dim, nodes.size()/dim);
  dumper.setPrefix("./");
  dumper.setVTUSubDirectory("test_tetra4-VTUs");
  dumper.init();
  dumper.dump("test_tetra4");
}
