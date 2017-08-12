/**
 * @file   test_data_distribution.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Sep 05 2014
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  Test the mesh distribution on creation of a distributed synchonizer
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
#include "aka_iterators.hh"
#include "aka_random_generator.hh"
#include "element_group.hh"
#include "element_synchronizer.hh"
#include "mesh_iterators.hh"
#include "mesh_partition_mesh_data.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize(argc, argv);

  const UInt spatial_dimension = 3;

  Mesh mesh_group_after(spatial_dimension, "after");
  Mesh mesh_group_before(spatial_dimension, "before");

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  if (prank == 0) {
    mesh_group_before.read("data_split.msh");
    mesh_group_after.read("data_split.msh");

    mesh_group_before.registerData<UInt>("global_id");
    mesh_group_after.registerData<UInt>("global_id");

    for (const auto & type : mesh_group_after.elementTypes(_all_dimensions)) {
      auto & gidb = mesh_group_before.getDataPointer<UInt>("global_id", type);
      auto & gida = mesh_group_after.getDataPointer<UInt>("global_id", type);

      for (auto && data : zip(arange(gida.getSize()), gida, gidb)) {
        std::get<1>(data) = std::get<0>(data);
        std::get<2>(data) = std::get<0>(data);
      }
    }
  }

  RandomGenerator<UInt>::seed(1);
  mesh_group_before.distribute();

  RandomGenerator<UInt>::seed(1);
  mesh_group_after.distribute();

  if (prank == 0)
    std::cout << mesh_group_after;

  for (const auto & agrp : ElementGroupsIterable(mesh_group_after)) {
    const ElementGroup & bgrp =
        mesh_group_before.getElementGroup(agrp.getName());

    for (auto && ghost_type : ghost_types) {
      for (const auto & type : bgrp.elementTypes(_all_dimensions, ghost_type)) {
        Array<UInt> & gidb = mesh_group_before.getDataPointer<UInt>(
            "global_id", type, ghost_type);
        Array<UInt> & gida = mesh_group_after.getDataPointer<UInt>(
            "global_id", type, ghost_type);

        Array<UInt> bgelem(bgrp.getElements(type, ghost_type));
        Array<UInt> agelem(agrp.getElements(type, ghost_type));

        std::transform(agelem.begin(), agelem.end(), agelem.begin(),
                       [&gida](auto & i) { return gida(i); });
        std::transform(bgelem.begin(), bgelem.end(), bgelem.begin(),
                       [&gidb](auto & i) { return gidb(i); });

        std::sort(bgelem.begin(), bgelem.end());
        std::sort(agelem.begin(), agelem.end());

        if (not std::equal(bgelem.begin(), bgelem.end(), agelem.begin())) {
          std::cerr << "The filters array for the group " << agrp.getName()
                    << " and for the element type " << type << ", "
                    << ghost_type << " do not match" << std::endl;

          debug::setDebugLevel(dblTest);

          std::cerr << bgelem << std::endl;
          std::cerr << agelem << std::endl;
          debug::debugger.exit(EXIT_FAILURE);
        }
      }
    }
  }

  for (const auto & agrp : NodeGroupsIterable(mesh_group_after)) {
    const NodeGroup & bgrp = mesh_group_before.getNodeGroup(agrp.getName());

    const Array<UInt> & gidb = mesh_group_before.getGlobalNodesIds();
    const Array<UInt> & gida = mesh_group_after.getGlobalNodesIds();

    Array<UInt> bgnode(0, 1);
    Array<UInt> agnode(0, 1);

    for (auto && pair : zip(bgrp, agrp)) {
      UInt a,b;
      std::tie(b, a) = pair;
      if (psize > 1) {
        if (mesh_group_before.isLocalOrMasterNode(b))
          bgnode.push_back(gidb(b));

        if (mesh_group_after.isLocalOrMasterNode(a))
          agnode.push_back(gida(a));
      }
    }

    std::sort(bgnode.begin(), bgnode.end());
    std::sort(agnode.begin(), agnode.end());

    if (!std::equal(bgnode.begin(), bgnode.end(), agnode.begin())) {
      std::cerr << "The filters array for the group " << agrp.getName() << " do not match"
                << std::endl;

      debug::setDebugLevel(dblTest);

      std::cerr << bgnode << std::endl;
      std::cerr << agnode << std::endl;
      debug::debugger.exit(EXIT_FAILURE);
    }
  }

  mesh_group_after.getElementGroup("inside").setBaseName("after_inside");
  mesh_group_after.getElementGroup("inside").dump();
  mesh_group_after.getElementGroup("outside").setBaseName("after_outside");
  mesh_group_after.getElementGroup("outside").dump();
  mesh_group_after.getElementGroup("volume").setBaseName("after_volume");
  mesh_group_after.getElementGroup("volume").dump();

  mesh_group_before.getElementGroup("inside").setBaseName("before_inside");
  mesh_group_before.getElementGroup("inside").dump();
  mesh_group_before.getElementGroup("outside").setBaseName("before_outside");
  mesh_group_before.getElementGroup("outside").dump();
  mesh_group_before.getElementGroup("volume").setBaseName("before_volume");
  mesh_group_before.getElementGroup("volume").dump();

  finalize();

  return EXIT_SUCCESS;
}
