/**
 * @file   dumper_low_level.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Aug 17 2015
 *
 * @brief  Example of dumper::DumperIOHelper low-level methods.
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "element_group.hh"
#include "group_manager.hh"
#include "mesh.hh"

#include "dumper_elemental_field.hh"
#include "dumper_nodal_field.hh"

#include "dumper_iohelper_paraview.hh"
#include "locomotive_tools.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char * argv[]) {

  /* This example aims at illustrating how to manipulate low-level methods of
     DumperIOHelper.
     The aims is to visualize a colorized moving train with Paraview */

  initialize(argc, argv);

  // To start let us load the swiss train mesh and its mesh data information.
  // We aknowledge here a weel-known swiss industry for mesh donation.
  UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("swiss_train.msh");

  Array<Real> & nodes = mesh.getNodes();
  UInt nb_nodes = mesh.getNbNodes();

  /* swiss_train.msh has the following physical groups that can be viewed with
    GMSH:
    "$MeshFormat
     2.2 0 8
     $EndMeshFormat
     $PhysicalNames
     6
     2 1 "red"
     2 2 "white"
     2 3 "lwheel_1"
     2 4 "lwheel_2"
     2 5 "rwheel_2"
     2 6 "rwheel_1"
     $EndPhysicalNames
     ..."
   */

  // Grouping nodes and elements belonging to train wheels (=four mesh data)
  ElementGroup & wheels_elements =
      mesh.createElementGroup("wheels", spatial_dimension);
  wheels_elements.append(mesh.getElementGroup("lwheel_1"));
  wheels_elements.append(mesh.getElementGroup("lwheel_2"));
  wheels_elements.append(mesh.getElementGroup("rwheel_1"));
  wheels_elements.append(mesh.getElementGroup("rwheel_2"));

  const Array<UInt> & lnode_1 =
      (mesh.getElementGroup("lwheel_1")).getNodeGroup().getNodes();
  const Array<UInt> & lnode_2 =
      (mesh.getElementGroup("lwheel_2")).getNodeGroup().getNodes();
  const Array<UInt> & rnode_1 =
      (mesh.getElementGroup("rwheel_1")).getNodeGroup().getNodes();
  const Array<UInt> & rnode_2 =
      (mesh.getElementGroup("rwheel_2")).getNodeGroup().getNodes();

  /* Note this Array is constructed with three components in order to warp train
     deformation on Paraview. A more appropriate way to do this is to set a
     padding in the NodalField (See example_dumpable_interface.cc.) */
  Array<Real> displacement(nb_nodes, 3);

  // ElementalField are constructed with an ElementTypeMapArray.
  ElementTypeMapArray<UInt> colour;
  colour.initialize(mesh, _with_nb_element = true);

  /* ------------------------------------------------------------------------ */
  /* Dumper creation                                                          */
  /* ------------------------------------------------------------------------ */

  // Creation of two DumperParaview. One for full mesh, one for a filtered
  // mesh.
  DumperParaview dumper("train", "./paraview/dumper", false);
  DumperParaview wheels("wheels", "./paraview/dumper", false);

  // Register the full mesh
  dumper.registerMesh(mesh);

  // Register a filtered mesh limited to nodes and elements from wheels groups
  wheels.registerFilteredMesh(mesh, wheels_elements.getElements(),
                              wheels_elements.getNodeGroup().getNodes());

  // Generate an output file of the two mesh registered.
  dumper.dump();
  wheels.dump();

  /* At this stage no fields are attached to the two dumpers.  To do so, a
     dumper::Field object has to be created.  Several types of dumper::Field
     exist. In this example we present two of them.

     NodalField to describe nodal displacements of our train.
     ElementalField handling the color of our different part.
  */

  // NodalField are constructed with an Array.
  auto displ_field = std::make_shared<dumper::NodalField<Real>>(displacement);
  auto colour_field = std::make_shared<dumper::ElementalField<UInt>>(colour);

  // Register the freshly created fields to our dumper.
  dumper.registerField("displacement", displ_field);
  dumper.registerField("colour", colour_field);

  // For the dumper wheels, fields have to be filtered at registration.
  // Filtered NodalField can be simply registered by adding an Array<UInt>
  // listing the nodes.
  auto displ_field_wheel = std::make_shared<dumper::NodalField<Real, true>>(
      displacement, 0, 0, &(wheels_elements.getNodeGroup().getNodes()));
  wheels.registerField("displacement", displ_field_wheel);

  // For the ElementalField, an ElementTypeMapArrayFilter has to be created.
  ElementTypeMapArrayFilter<UInt> filtered_colour(
      colour, wheels_elements.getElements());

  auto colour_field_wheel =
      std::make_shared<dumper::ElementalField<UInt, Vector, true>>(
          filtered_colour);
  wheels.registerField("colour", colour_field_wheel);

  /* ------------------------------------------------------------------------ */
  // Now that the dumpers are created and the fields are associated, let's
  // paint and move the train!

  // Fill the ElementTypeMapArray colour according to mesh data information.
  fillColour(mesh, colour);

  // Apply displacement and wheels rotation.
  Real tot_displacement = 50.;
  Real radius = 1.;
  UInt nb_steps = 100;
  Real theta = tot_displacement / radius;

  Vector<Real> l_center(3);
  Vector<Real> r_center(3);

  for (UInt i = 0; i < spatial_dimension; ++i) {
    l_center(i) = nodes(14, i);
    r_center(i) = nodes(2, i);
  }

  for (UInt i = 0; i < nb_steps; ++i) {
    displacement.clear();

    Real angle = (Real)i / (Real)nb_steps * theta;
    applyRotation(l_center, angle, nodes, displacement, lnode_1);
    applyRotation(l_center, angle, nodes, displacement, lnode_2);
    applyRotation(r_center, angle, nodes, displacement, rnode_1);
    applyRotation(r_center, angle, nodes, displacement, rnode_2);

    for (UInt j = 0; j < nb_nodes; ++j) {
      displacement(j, 0) += (Real)i / (Real)nb_steps * tot_displacement;
    }
    // Output results after each moving steps for main and wheel dumpers.
    dumper.dump();
    wheels.dump();
  }

  finalize();

  return 0;
}
