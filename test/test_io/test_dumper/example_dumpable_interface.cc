/**
 * @file   example_dumpable_interface.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Jul  2 14:34:41 2015
 *
 * @brief  Example of IOHelper dumper advanced usage  
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
#include "element_group.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "dumper_paraview.hh"
#include "locomotive.hh" 

#define PI 3.141592653589
/* -------------------------------------------------------------------------- */
using namespace akantu;


int main(int argc, char *argv[]) {

  /// In this example, we 
  initialize(argc, argv);

  /// To start let us load the swiss train mesh with its mesh data
  UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("swiss_train.msh");    
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  /* swiss_train.msh has the following physical groups that can be viewed with GMSH:
    "$MeshFormat
     2.2 0 8
     $EndMeshFormat
     $PhysicalNames
     6
     2 1 "rouge"
     2 2 "blanc"
     2 3 "lwheel_1"
     2 4 "lwheel_2"
     2 5 "rwheel_2"
     2 6 "rwheel_1"
     $EndPhysicalNames
     ..."
   */

  /// Create grouping nodes and elements belonging to train wheels (=four mesh data)
  NodeGroup wheels_nodes("wheels_nodes", mesh);
  wheels_nodes.append(mesh.getElementGroup("lwheel_1").getNodeGroup());
  wheels_nodes.append(mesh.getElementGroup("lwheel_2").getNodeGroup());
  wheels_nodes.append(mesh.getElementGroup("rwheel_1").getNodeGroup());
  wheels_nodes.append(mesh.getElementGroup("rwheel_2").getNodeGroup());
 
  ElementGroup & wheels_elements = mesh.createElementGroup("wheels_elements", spatial_dimension, wheels_nodes);
  wheels_elements.append(mesh.getElementGroup("lwheel_1"));
  wheels_elements.append(mesh.getElementGroup("lwheel_2"));
  wheels_elements.append(mesh.getElementGroup("rwheel_1"));
  wheels_elements.append(mesh.getElementGroup("rwheel_2"));

  DumperParaview dumper("train", "paraview/dumpable", false);
  mesh.registerExternalDumper(dumper, "train", true);
  mesh.addDumpMesh(mesh);
  
  mesh.registerExternalDumper(mesh.getGroupDumper("paraview_wheels_elements", "wheels_elements"),"wheels_elements");
  mesh.setDirectoryToDumper("wheels_elements","paraview/dumpable");

  Array<Real> & node = mesh.getNodes();
  UInt nb_nodes = mesh.getNbNodes();
  Array<Real> * displacement = new Array<Real>(nb_nodes,spatial_dimension);
  UInt padding_size = 3;

  ElementTypeMapArray<UInt> & colour = mesh.registerData<UInt>("colour");
  mesh.initElementTypeMapArray(colour, 1, spatial_dimension, false, _ek_regular, true);   
  
  mesh.addDumpFieldExternal("displacement", mesh.createNodalField(displacement, "all",padding_size));
  mesh.addDumpFieldExternal("color", mesh.createFieldFromAttachedData<UInt>("colour", "all", _ek_regular)); 

  mesh.addDumpFieldExternalToDumper("wheels_elements", "displacement", 
				    mesh.createNodalField(displacement, "wheels_elements",padding_size));
  mesh.addDumpFieldExternalToDumper("wheels_elements", "colour", 
				    mesh.createFieldFromAttachedData<UInt>("colour", "wheels_elements", _ek_regular));  
  
  /// Fill the ElementTypeMapArray colour according to mesh data information
  ElementTypeMapArray<std::string> * phys_data = &(mesh.getData<std::string>("physical_names"));
  Array<std::string> & txt_colour = (*phys_data)(_triangle_3);
  Array<UInt> & id_colour = (colour)(_triangle_3);

  for (UInt i = 0; i < txt_colour.getSize(); ++i) {
  
    std::string phy_name = txt_colour[i];

    if (phy_name == "rouge")
      id_colour[i] = 3;
    else if (phy_name == "blanc"||phy_name == "lwheel_1"||phy_name == "rwheel_1")
      id_colour[i] = 2;
    else
      id_colour[i] = 1;
  }
 
  /// Apply displacement and wheels rotation
  Real tot_displacement = 50.;
  Real radius = 1.;
  UInt nb_steps = 500;
  Real theta = tot_displacement/radius;

  Array<Real> l_center(spatial_dimension);
  Array<Real> r_center(spatial_dimension);

  const Array<UInt> & lnode_1 = (mesh.getElementGroup("lwheel_1")).getNodes();  
  const Array<UInt> & lnode_2 = (mesh.getElementGroup("lwheel_2")).getNodes();  
  const Array<UInt> & rnode_1 = (mesh.getElementGroup("rwheel_1")).getNodes();  
  const Array<UInt> & rnode_2 = (mesh.getElementGroup("rwheel_2")).getNodes();  

  for (UInt i = 0; i < spatial_dimension; ++i) {

    l_center(i) = node(14,i);
    r_center(i) = node(2,i);
  }

  for (UInt i = 0; i < nb_steps; ++i) {
    
    displacement->clear();
  
    Real angle = (Real)i / (Real)nb_steps * theta;
    applyRotation(l_center, angle, node, displacement, lnode_1);
    applyRotation(l_center, angle, node, displacement, lnode_2);
    applyRotation(r_center, angle, node, displacement, rnode_1);
    applyRotation(r_center, angle, node, displacement, rnode_2);

    for (UInt j = 0; j < nb_nodes; ++j) {
      (*displacement)(j,0) += (Real)i / (Real)nb_steps *tot_displacement;
    }
    /// Results are output after each moving steps for main and wheel dumpers.  
    mesh.dump();
    mesh.dump("wheels_elements");
  }
  return 0;
}
