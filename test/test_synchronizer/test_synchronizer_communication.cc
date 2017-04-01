/**
 * @file   test_synchronizer_communication.cc
 *
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 01 2010
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  test to synchronize barycenters
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "aka_common.hh"
#include "aka_random_generator.hh"
#include "element_synchronizer.hh"
#include "mesh.hh"
#include "mesh_partition_scotch.hh"
#include "synchronizer_registry.hh"
/* -------------------------------------------------------------------------- */

//#define DUMP

#if defined(AKANTU_USE_IOHELPER) && defined(DUMP)
#include "dumper_paraview.hh"
#endif // AKANTU_USE_IOHELPER

#include "test_data_accessor.hh"

using namespace akantu;

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);

  UInt spatial_dimension = 3;

  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int prank = comm.whoAmI();
  bool wait = true;
  if (argc > 1) {
    if (prank == 0)
      while (wait)
        ;
  }

  if (prank == 0)
    mesh.read("cube.msh");

  mesh.distribute();

  /// compute barycenter for each facet
  ElementTypeMapArray<Real> barycenters("barycenters", "", 0);
  mesh.initElementTypeMapArray(barycenters, spatial_dimension,
                               spatial_dimension);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last_type =
        mesh.lastType(spatial_dimension, ghost_type);

    for (; it != last_type; ++it) {
      UInt nb_element = mesh.getNbElement(*it, ghost_type);
      Array<Real> & barycenter = barycenters(*it, ghost_type);
      barycenter.resize(nb_element);

      Array<Real>::iterator<Vector<Real>> bary_it =
          barycenter.begin(spatial_dimension);

      for (UInt elem = 0; elem < nb_element; ++elem, ++bary_it)
        mesh.getBarycenter(elem, *it, bary_it->storage(), ghost_type);
    }
  }

  AKANTU_DEBUG_INFO("Creating TestAccessor");
  TestAccessor test_accessor(mesh, barycenters);
  SynchronizerRegistry synch_registry;
  synch_registry.registerDataAccessor(test_accessor);
  synch_registry.registerSynchronizer(mesh.getElementSynchronizer(), _gst_test);

  AKANTU_DEBUG_INFO("Synchronizing tag");
  synch_registry.synchronize(_gst_test);

  // Checking the tags in MeshData (not a very good test because they're all
  // identical,
  // but still...)
  Mesh::type_iterator it = mesh.firstType(_all_dimensions);
  Mesh::type_iterator last_type = mesh.lastType(_all_dimensions);

  for (; it != last_type; ++it) {
    Array<UInt> & tags = mesh.getData<UInt>("tag_0", *it);
    Array<UInt>::const_vector_iterator tags_it = tags.begin(1);
    Array<UInt>::const_vector_iterator tags_end = tags.end(1);
    AKANTU_DEBUG_ASSERT(
        mesh.getNbElement(*it) == tags.getSize(),
        "The number of tags does not match the number of elements on rank "
            << prank << ".");
    std::cout << std::dec << " I am rank " << prank
              << " and here's my MeshData dump for types " << *it
              << " (it should contain " << mesh.getNbElement(*it)
              << " elements and it has " << tags.getSize()
              << "!) :" << std::endl;
    std::cout << std::hex;

    debug::setDebugLevel(dblTest);
    for (; tags_it != tags_end; ++tags_it) {
      std::cout << tags_it->operator()(0) << " ";
    }

    debug::setDebugLevel(dblInfo);
    std::cout << std::endl;
  }

#if defined(AKANTU_USE_IOHELPER) && defined(DUMP)
  DumperParaview dumper("test-scotch-partition");
  dumper.registerMesh(mesh, spatial_dimension, _not_ghost);
  dumper.registerField("partitions",
   		       new DumperIOHelper::ElementPartitionField<>(mesh,
                                                                   spatial_dimension, _not_ghost));
  dumper.dump();

  DumperParaview dumper_ghost("test-scotch-partition-ghost");
  dumper_ghost.registerMesh(mesh, spatial_dimension, _ghost);
  dumper_ghost.registerField("partitions",
   			     new DumperIOHelper::ElementPartitionField<>(mesh,
                                                                         spatial_dimension, _ghost));
  dumper_ghost.dump();
#endif //AKANTU_USE_IOHELPER

  finalize();

  return EXIT_SUCCESS;
}
