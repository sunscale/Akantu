/**
 * @file   test_mesh_partitionate_scotch.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Sep 12 2010
 * @date last modification: Mon Jan 22 2018
 *
 * @brief  test of internal facet extraction
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_partition_scotch.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#include "dumper_elemental_field.hh"
#include "dumper_iohelper_paraview.hh"
#endif // AKANTU_USE_IOHELPER

using namespace akantu;
/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);
  debug::setDebugLevel(akantu::dblDump);

  int dim = 2;

  Mesh mesh(dim);
  mesh.read("triangle.msh");

  MeshPartitionScotch partition(mesh, dim);
  partition.partitionate(8);

#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper("test-scotch-partition");
  auto field = std::make_shared<dumpers::ElementalField<UInt>>(
      partition.getPartitions(), dim);
  dumper.registerMesh(mesh, dim);
  dumper.registerField("partitions", field);
  dumper.dump();
#endif // AKANTU_USE_IOHELPER

  partition.reorder();
  mesh.write("triangle_reorder.msh");

  finalize();

  return EXIT_SUCCESS;
}
