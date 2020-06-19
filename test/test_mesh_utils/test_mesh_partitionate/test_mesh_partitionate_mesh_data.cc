/**
 * @file   test_mesh_partitionate_mesh_data.cc
 *
 * @author Dana Christen <dana.christen@epfl.ch>
 *
 * @date creation: Wed May 08 2013
 * @date last modification: Tue Nov 07 2017
 *
 * @brief  test of manual partitioner
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_partition_mesh_data.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#include "dumper_elemental_field.hh"
#include "dumper_paraview.hh"
#endif // AKANTU_USE_IOHELPER
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);

  UInt dim = 2;
  UInt nb_partitions = 8;
  akantu::Mesh mesh(dim);
  mesh.read("quad.msh");

  ElementTypeMapArray<UInt> partition;
  UInt nb_component = 1;

  GhostType gt = _not_ghost;
  for (auto & type : mesh.elementTypes(dim, gt)) {
    UInt nb_element = mesh.getNbElement(type, gt);
    partition.alloc(nb_element, nb_component, type, gt);
    Array<UInt> & type_partition_reference = partition(type, gt);
    for (UInt i(0); i < nb_element; ++i) {
      Vector<Real> barycenter(dim);
      Element element{type, i, gt};
      mesh.getBarycenter(element, barycenter);

      Real real_proc = barycenter[0] * nb_partitions;
      if (std::abs(real_proc - round(real_proc)) <
          10 * std::numeric_limits<Real>::epsilon()) {
        type_partition_reference(i) = round(real_proc);
      } else {
        std::cout << "*";
        type_partition_reference(i) = floor(real_proc);
      }
      std::cout << "Assigned proc " << type_partition_reference(i)
                << " to elem " << i << " (type " << type
                << ", barycenter x-coordinate " << barycenter[0] << ")"
                << std::endl;
    }
  }

  akantu::MeshPartitionMeshData * partitioner =
      new akantu::MeshPartitionMeshData(mesh, dim);
  partitioner->setPartitionMapping(partition);
  partitioner->partitionate(nb_partitions);

  for (auto & type : mesh.elementTypes(dim, gt)) {
    UInt nb_element = mesh.getNbElement(type, gt);
    const Array<UInt> & type_partition_reference = partition(type, gt);
    const Array<UInt> & type_partition = partitioner->getPartitions()(type, gt);
    for (UInt i(0); i < nb_element; ++i) {
      if (not(type_partition(i) == type_partition_reference(i))) {
        std::cout << "Incorrect partitioning" << std::endl;
        return 1;
      }
    }
  }

#ifdef DEBUG_TEST
  DumperParaview dumper("test-mesh-data-partition");
  dumpers::Field * field1 =
      new dumpers::ElementalField<UInt>(partitioner->getPartitions(), dim);
  dumpers::Field * field2 = new dumpers::ElementalField<UInt>(partition, dim);
  dumper.registerMesh(mesh, dim);
  dumper.registerField("partitions", field1);
  dumper.registerField("partitions_ref", field2);
  dumper.dump();
#endif

  delete partitioner;

  finalize();

  return EXIT_SUCCESS;
}
