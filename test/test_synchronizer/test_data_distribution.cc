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
#include "aka_common.hh"
#include "element_group.hh"
#include "element_synchronizer.hh"
#include "mesh_partition_mesh_data.hh"
#include "aka_random_generator.hh"

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

    for (Mesh::type_iterator tit = mesh_group_after.firstType(_all_dimensions);
         tit != mesh_group_after.lastType(_all_dimensions); ++tit) {
      Array<UInt> & gidb =
          mesh_group_before.getDataPointer<UInt>("global_id", *tit);
      Array<UInt> & gida =
          mesh_group_after.getDataPointer<UInt>("global_id", *tit);

      Array<UInt>::scalar_iterator ait = gida.begin();
      Array<UInt>::scalar_iterator bit = gidb.begin();
      Array<UInt>::scalar_iterator end = gida.end();
      for (UInt i = 0; ait != end; ++ait, ++i, ++bit) {
        *bit = i;
        *ait = i;
      }
    }
  }

  RandomGenerator<UInt>::seed(1);
  mesh_group_before.distribute();

  RandomGenerator<UInt>::seed(1);
  mesh_group_after.distribute();

  if (prank == 0)
    std::cout << mesh_group_after;

  GroupManager::element_group_iterator grp_ait =
      mesh_group_after.element_group_begin();
  GroupManager::element_group_iterator grp_end =
      mesh_group_after.element_group_end();
  for (; grp_ait != grp_end; ++grp_ait) {
    std::string grp = grp_ait->first;
    const ElementGroup & bgrp = mesh_group_before.getElementGroup(grp);
    const ElementGroup & agrp = *grp_ait->second;

    for (ghost_type_t::iterator git = ghost_type_t::begin();
         git != ghost_type_t::end(); ++git) {
      GhostType ghost_type = *git;

      for (Mesh::type_iterator tit =
               bgrp.firstType(_all_dimensions, ghost_type);
           tit != bgrp.lastType(_all_dimensions, ghost_type); ++tit) {
        Array<UInt> & gidb = mesh_group_before.getDataPointer<UInt>(
            "global_id", *tit, ghost_type);
        Array<UInt> & gida = mesh_group_after.getDataPointer<UInt>(
            "global_id", *tit, ghost_type);

        Array<UInt> bgelem(bgrp.getElements(*tit, ghost_type));
        Array<UInt> agelem(agrp.getElements(*tit, ghost_type));

        Array<UInt>::scalar_iterator ait = agelem.begin();
        Array<UInt>::scalar_iterator bit = bgelem.begin();
        Array<UInt>::scalar_iterator end = agelem.end();
        for (; ait != end; ++ait, ++bit) {
          *bit = gidb(*bit);
          *ait = gida(*ait);
        }

        std::sort(bgelem.begin(), bgelem.end());
        std::sort(agelem.begin(), agelem.end());

        if (!std::equal(bgelem.begin(), bgelem.end(), agelem.begin())) {
          std::cerr << "The filters array for the group " << grp
                    << " and for the element type " << *tit << ", "
                    << ghost_type << " do not match" << std::endl;

          debug::setDebugLevel(dblTest);

          std::cerr << bgelem << std::endl;
          std::cerr << agelem << std::endl;
          debug::debugger.exit(EXIT_FAILURE);
        }
      }
    }
  }

  GroupManager::node_group_iterator ngrp_ait =
      mesh_group_after.node_group_begin();
  GroupManager::node_group_iterator ngrp_end =
      mesh_group_after.node_group_end();
  for (; ngrp_ait != ngrp_end; ++ngrp_ait) {
    std::string grp = ngrp_ait->first;
    const NodeGroup & bgrp = mesh_group_before.getNodeGroup(grp);
    const NodeGroup & agrp = *ngrp_ait->second;

    const Array<UInt> & gidb = mesh_group_before.getGlobalNodesIds();
    const Array<UInt> & gida = mesh_group_after.getGlobalNodesIds();

    Array<UInt> bgnode(0, 1);
    Array<UInt> agnode(0, 1);

    Array<UInt>::const_scalar_iterator ait = agrp.begin();
    Array<UInt>::const_scalar_iterator bit = bgrp.begin();
    Array<UInt>::const_scalar_iterator end = agrp.end();
    for (; ait != end; ++ait, ++bit) {
      if (psize > 1) {
        if (mesh_group_before.isLocalOrMasterNode(*bit))
          bgnode.push_back(gidb(*bit));

        if (mesh_group_after.isLocalOrMasterNode(*ait))
          agnode.push_back(gida(*ait));
      }
    }

    std::sort(bgnode.begin(), bgnode.end());
    std::sort(agnode.begin(), agnode.end());

    if (!std::equal(bgnode.begin(), bgnode.end(), agnode.begin())) {
      std::cerr << "The filters array for the group " << grp << " do not match"
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
