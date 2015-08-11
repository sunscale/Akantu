/**
 * @file   swiss_train.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Jul  2 14:34:41 2015
 *
 * @brief  Example of IOHelper dumper
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
#include "dumper_elemental_field.hh"
#include "dumper_nodal_field.hh" 

#define PI 3.141592653589
/* -------------------------------------------------------------------------- */
using namespace akantu;

/* -------------------------------------------------------------------------- */
void applyRotation(Array<Real> center, Real angle, Array<Real> & nodes,
		   Array<Real> * displacement, Array<UInt> node_group) {

  for (UInt i = 0; i < node_group.getSize(); ++i) {
    
    Array<Real> pos_rel(center.getSize());

    for (UInt j = 0; j < center.getSize(); ++j) {
      
      pos_rel(j) = nodes(node_group(i),j) - center(j);
    }

    Real radius = std::sqrt(pos_rel[0]*pos_rel[0]+pos_rel[1]*pos_rel[1]);

    if (std::abs(radius) < Math::getTolerance()) continue;

    Real phi_i = std::acos(pos_rel[0]/radius);

    if (pos_rel[1] < 0) phi_i *= -1;

    (*displacement)(node_group(i),0) = std::cos(phi_i-angle)*radius - pos_rel[0];
    (*displacement)(node_group(i),1) = std::sin(phi_i-angle)*radius - pos_rel[1];
  }
}

int main(int argc, char *argv[]) {

  initialize(argc, argv);
 
  UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("swiss_train.msh");    
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  Real tot_displacement = 50.;
  Real radius = 1.;
  UInt nb_steps = 500;
  Real theta = tot_displacement/radius;

  DumperParaview dumper("barras_", "./paraview", false);
  DumperParaview wheels("barras_", "./paraview", false);

  dumper.registerMesh(mesh);

  dumper.setDirectory("paraview/train/");
  wheels.setDirectory("paraview/train/");
  dumper.setBaseName("swiss_train");
  wheels.setBaseName("wheels");

  //mesh.registerExternalDumper(dumper, "mydumper");
  dumper.dump();

  Array<Real> & node = mesh.getNodes();
  UInt nb_nodes = mesh.getNbNodes();

  Array<Real> * displacement = new Array<Real>(nb_nodes,3);
  
  ElementTypeMapArray<std::string> * phys_data = &(mesh.getData<std::string>("physical_names"));
  ElementTypeMapArray<UInt> colour;
  mesh.initElementTypeMapArray(colour, 1, spatial_dimension, false, _ek_regular, true); 
  //ETMArray to initialize, nb_component, spatial_dim, nb_node_by_elem_size?, type, allocate size getNbElement(type)?

  dumper::Field * field = new dumper::NodalField<Real>(*displacement);
  dumper::Field * color = new dumper::ElementalField<UInt>(colour);

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


  dumper.registerField("displacement", field);
  dumper.registerField("colour", color);
  
  const Array<UInt> & lnode_1 = (mesh.getElementGroup("lwheel_1")).getNodes();  
  const Array<UInt> & lnode_2 = (mesh.getElementGroup("lwheel_2")).getNodes();  
  const Array<UInt> & rnode_1 = (mesh.getElementGroup("rwheel_1")).getNodes();  
  const Array<UInt> & rnode_2 = (mesh.getElementGroup("rwheel_2")).getNodes();  
  
  NodeGroup wheels_nodes("wheels_nodes", mesh);
  wheels_nodes.append(mesh.getElementGroup("lwheel_1").getNodeGroup());
  wheels_nodes.append(mesh.getElementGroup("lwheel_2").getNodeGroup());
  wheels_nodes.append(mesh.getElementGroup("rwheel_1").getNodeGroup());
  wheels_nodes.append(mesh.getElementGroup("rwheel_2").getNodeGroup());
  ElementGroup wheels_elements("wheels_elements", mesh, wheels_nodes);
  wheels_elements.append(mesh.getElementGroup("lwheel_1"));
  wheels_elements.append(mesh.getElementGroup("lwheel_2"));
  wheels_elements.append(mesh.getElementGroup("rwheel_1"));
  wheels_elements.append(mesh.getElementGroup("rwheel_2"));

  wheels.registerFilteredMesh(mesh,
			      wheels_elements.getElements(),
			      wheels_elements.getNodes());

  wheels.registerField("displacement", 
		       new dumper::NodalField<Real,true>(*displacement,0,0,&(wheels_elements.getNodes())));

  const ElementTypeMapArray<UInt> & elemental_filter = wheels_elements.getElements();
  ElementTypeMapArrayFilter<UInt> * filter = new ElementTypeMapArrayFilter<UInt>(colour,elemental_filter);

  wheels.registerField("colour", 
		       new dumper::ElementalField<UInt,Vector,true>(*filter));


  Array<Real> l_center(spatial_dimension);
  Array<Real> r_center(spatial_dimension);

  for (UInt i = 0; i < spatial_dimension; ++i) {

    l_center(i) = node(14,i);
    r_center(i) = node(2,i);
  }

  std::cout << l_center[0] << " " << l_center[1] << std::endl;
  std::cout << r_center[0] << " " << r_center[1] << std::endl;

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
    dumper.dump();
    wheels.dump();
  }
  return 0;
}
