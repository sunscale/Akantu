/**
 * @file   example_dumper_low_level.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Jul  2 14:34:41 2015
 *
 * @brief  Example of dumper::DumperIOHelper low-level methods.   
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
#include "locomotive.hh" 

#define PI 3.141592653589
/* -------------------------------------------------------------------------- */
using namespace akantu;


int main(int argc, char *argv[]) {

  /* This example aims at illustrating how to manipulate low-level methods of DumperIOHelper.
     The aims is to visualize a colorized moving train with Paraview */

  initialize(argc, argv);

  /// To start let us load the swiss train mesh and its mesh data information.
  /// We aknowledge here a weel-known swiss industry for mesh donation.
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

  /// Creation of two DumperParaview. One for full mesh, one for a filtered mesh.
  DumperParaview dumper("train", "./paraview/dumper", false);
  DumperParaview wheels("wheels", "./paraview/dumper", false);

  /// Register the full mesh
  dumper.registerMesh(mesh);

  /// Grouping nodes and elements belonging to train wheels (=four mesh data)
  NodeGroup & wheels_nodes = mesh.createNodeGroup("wheels_nodes");
  wheels_nodes.append(mesh.getElementGroup("lwheel_1").getNodeGroup());
  wheels_nodes.append(mesh.getElementGroup("lwheel_2").getNodeGroup());
  wheels_nodes.append(mesh.getElementGroup("rwheel_1").getNodeGroup());
  wheels_nodes.append(mesh.getElementGroup("rwheel_2").getNodeGroup());
  
  ElementGroup & wheels_elements = mesh.createElementGroup("wheels_elements", spatial_dimension, wheels_nodes);
  wheels_elements.append(mesh.getElementGroup("lwheel_1"));
  wheels_elements.append(mesh.getElementGroup("lwheel_2"));
  wheels_elements.append(mesh.getElementGroup("rwheel_1"));
  wheels_elements.append(mesh.getElementGroup("rwheel_2"));

  /// Register a filtered mesh limited to nodes and elements from wheels groups
  wheels.registerFilteredMesh(mesh,
			      wheels_elements.getElements(),
			      wheels_elements.getNodes());

  /// Generate an output file of the two mesh registered.  
  dumper.dump();
  wheels.dump();

  /// At this stage no fields are attached to the two meshes. 
  /// To do so, a dumper::Field object has to be created.
  /// Several types of dumper::Field exist. In this example we present two of them.

  /// NodalField to describe nodal displacements of our train.
  /// ElementalField handling the color of our different part.

  Array<Real> & node = mesh.getNodes();
  UInt nb_nodes = mesh.getNbNodes();

  /// NodalField are constructed with an Array.
  /* Note that this Array is constructed with three components
     in order to warp train deformation on Paraview.
     A more appropriate way to do this is to set a padding in the NodalField
     (See example_dumpable_interface.cc.) */
  Array<Real> * displacement = new Array<Real>(nb_nodes,3);
  dumper::Field * displ = new dumper::NodalField<Real>(*displacement);

  /// ElementalField are constructed with an ElementTypeMapArray. 
  ElementTypeMapArray<UInt> colour;
  mesh.initElementTypeMapArray(colour, 1, spatial_dimension, false, _ek_regular, true);   
  dumper::Field * color = new dumper::ElementalField<UInt>(colour);

  /// Register the freshly created fields to our dumper.
  dumper.registerField("displacement", displ);
  dumper.registerField("colour", color);

  /// For the dumper wheels, fields have to be filtered at registration.

  /// Filter NodalField can be simply registered by adding an Array<UInt> listing filtered nodes.    
  wheels.registerField("displacement", 
		       new dumper::NodalField<Real,true>(*displacement,0,0,&(wheels_elements.getNodes())));

  /// For ElementalField, an ElementTypeMapArrayFilter has to be created.
  const ElementTypeMapArray<UInt> & elemental_filter = wheels_elements.getElements();
  ElementTypeMapArrayFilter<UInt> * filter = new ElementTypeMapArrayFilter<UInt>(colour,elemental_filter);

  wheels.registerField("colour", 
		       new dumper::ElementalField<UInt,Vector,true>(*filter));

  ///Now that our dumpers are created and the two fields are associated, let's paint and move the train!

  /// Fill the ElementTypeMapArray colour according to mesh data information.
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

  /// Apply displacement and wheels rotation.
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
    /// Output results after each moving steps for main and wheel dumpers.  
    dumper.dump();
    wheels.dump();
  }
  return 0;
}
