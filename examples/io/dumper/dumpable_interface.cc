/**
 * @file   dumpable_interface.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Aug 17 2015
 * @date last modification: Mon Aug 31 2015
 *
 * @brief  Example of dumper::Dumpable interface.
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
#include "group_manager_inline_impl.cc"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include "dumpable_inline_impl.hh"
#include "dumper_iohelper_paraview.hh"
/* -------------------------------------------------------------------------- */
#include "locomotive_tools.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char * argv[]) {

  /*
    In this example, we present dumper::Dumpable which is an interface
    for other classes who want to dump themselves.
    Several classes of Akantu inheritate from Dumpable (Model, Mesh, ...).
    In this example we reproduce the same tasks as example_dumper_low_level.cc
    using this time Dumpable interface inherted by Mesh, NodeGroup and
    ElementGroup.
    It is then advised to read first example_dumper_low_level.cc.
  */

  initialize(argc, argv);

  // To start let us load the swiss train mesh and its mesh data information.
  UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("swiss_train.msh");

  /*
    swiss_train.msh has the following physical groups that can be viewed with
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

  // Grouping nodes and elements belonging to train wheels (=four mesh data).
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

  Array<Real> & node = mesh.getNodes();
  UInt nb_nodes = mesh.getNbNodes();

  // This time a 2D Array is created and a padding size of 3 is passed to
  // NodalField in order to warp train deformation on Paraview.
  Array<Real> displacement(nb_nodes, spatial_dimension);

  // Create an ElementTypeMapArray for the colour
  ElementTypeMapArray<UInt> colour("colour");
  colour.initialize(mesh, _with_nb_element = true);

  /* ------------------------------------------------------------------------ */
  /* Creating dumpers                                                         */
  /* ------------------------------------------------------------------------ */

  // Create dumper for the complete mesh and register it as default dumper.
  DumperParaview dumper("train", "./paraview/dumpable", false);
  mesh.registerExternalDumper(dumper, "train", true);
  mesh.addDumpMesh(mesh);

  // The dumper for the filtered mesh can be directly taken from the
  // ElementGroup and then registered as "wheels_elements" dumper.
  DumperIOHelper & wheels = mesh.getGroupDumper("paraview_wheels", "wheels");

  mesh.registerExternalDumper(wheels, "wheels");
  mesh.setDirectoryToDumper("wheels", "./paraview/dumpable");

  // Arrays and ElementTypeMapArrays can be added as external fields directly
  mesh.addDumpFieldExternal("displacement", displacement);

  ElementTypeMapArrayFilter<UInt> filtered_colour(
      colour, wheels_elements.getElements());

  auto colour_field_wheel =
      std::make_shared<dumper::ElementalField<UInt, Vector, true>>(
          filtered_colour);
  mesh.addDumpFieldExternal("color", colour_field_wheel);

  mesh.addDumpFieldExternalToDumper("wheels", "displacement", displacement);
  mesh.addDumpFieldExternalToDumper("wheels", "colour", colour);

  // For some specific cases the Fields should be created, as when you want to
  // pad an array
  auto displacement_vector_field =
      mesh.createNodalField(&displacement, "all", 3);
  mesh.addDumpFieldExternal("displacement_as_paraview_vector",
                            displacement_vector_field);
  mesh.addDumpFieldExternalToDumper("wheels", "displacement_as_paraview_vector",
                                    displacement_vector_field);

  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */

  // Fill the ElementTypeMapArray colour.
  fillColour(mesh, colour);

  /// Apply displacement and wheels rotation.
  Real tot_displacement = 50.;
  Real radius = 1.;
  UInt nb_steps = 100;
  Real theta = tot_displacement / radius;

  Vector<Real> l_center(spatial_dimension);
  Vector<Real> r_center(spatial_dimension);

  for (UInt i = 0; i < spatial_dimension; ++i) {
    l_center(i) = node(14, i);
    r_center(i) = node(2, i);
  }

  for (UInt i = 0; i < nb_steps; ++i) {
    displacement.clear();

    Real step_ratio = Real(i) / Real(nb_steps);
    Real angle = step_ratio * theta;

    applyRotation(l_center, angle, node, displacement, lnode_1);
    applyRotation(l_center, angle, node, displacement, lnode_2);
    applyRotation(r_center, angle, node, displacement, rnode_1);
    applyRotation(r_center, angle, node, displacement, rnode_2);

    for (UInt j = 0; j < nb_nodes; ++j) {
      displacement(j, _x) += step_ratio * tot_displacement;
    }
    /// Dump call is finally made through Dumpable interface.
    mesh.dump();
    mesh.dump("wheels");
  }

  finalize();

  return 0;
}
