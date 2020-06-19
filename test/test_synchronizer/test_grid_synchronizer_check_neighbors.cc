/**
 * @file   test_grid_synchronizer_check_neighbors.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Mar 11 2013
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Test the generation of neighbors list based on a akaentu::Grid
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

#include <fstream>
#include <iostream>
#include <string>

#include "aka_common.hh"
#include "communicator.hh"

using namespace akantu;

const UInt spatial_dimension = 3;

#include "test_grid_tools.hh"

void readNeighbors(
    std::ifstream & nin,
    neighbors_map_t<spatial_dimension>::type & neighbors_map_read) {
  std::string line;
  while (std::getline(nin, line)) {
    std::getline(nin, line);
    std::istringstream iss(line);
    UInt nb_neig;
    iss >> nb_neig;
    std::getline(nin, line);
    Point<spatial_dimension> pt;
    pt.read(line);
    std::getline(nin, line);
    for (UInt i = 0; i < nb_neig; ++i) {
      std::getline(nin, line);
      Point<spatial_dimension> ne;
      ne.read(line);
      neighbors_map_read[pt].push_back(ne);
    }
  }
}

int main(int argc, char * argv[]) {
  initialize(argc, argv);

  const auto & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  std::string file_ref = "neighbors_ref_1_0";
  std::string file = "neighbors_ref";
  std::stringstream sstr;
  sstr << file << "_" << psize << "_" << prank;
  file = sstr.str();

  std::ifstream nin;

  neighbors_map_t<spatial_dimension>::type neighbors_map_read;
  nin.open(file_ref.c_str());
  readNeighbors(nin, neighbors_map_read);
  nin.close();

  neighbors_map_t<spatial_dimension>::type neighbors_map;
  nin.open(file.c_str());
  readNeighbors(nin, neighbors_map);
  nin.close();

  neighbors_map_t<spatial_dimension>::type::iterator it_n =
      neighbors_map.begin();
  neighbors_map_t<spatial_dimension>::type::iterator end_n =
      neighbors_map.end();
  for (; it_n != end_n; ++it_n) {
    std::sort(it_n->second.begin(), it_n->second.end());

    std::vector<Point<spatial_dimension>>::iterator it_v = it_n->second.begin();
    std::vector<Point<spatial_dimension>>::iterator end_v = it_n->second.end();

    neighbors_map_t<spatial_dimension>::type::iterator it_nr =
        neighbors_map_read.find(it_n->first);
    if (it_nr == neighbors_map_read.end())
      AKANTU_ERROR(
          "Argh what is this point that is not present in the ref file "
          << it_n->first);

    std::vector<Point<spatial_dimension>>::iterator it_vr =
        it_nr->second.begin();
    std::vector<Point<spatial_dimension>>::iterator end_vr =
        it_nr->second.end();

    for (; it_v != end_v && it_vr != end_vr; ++it_v, ++it_vr) {
      if (*it_vr != *it_v)
        AKANTU_ERROR("Neighbors does not match " << *it_v << " != " << *it_vr
                                                 << " neighbor of "
                                                 << it_n->first);
    }

    if (it_v == end_v && it_vr != end_vr) {
      AKANTU_ERROR("Some neighbors of " << it_n->first << " are missing!");
    }

    if (it_v != end_v && it_vr == end_vr)
      AKANTU_ERROR("Some neighbors of " << it_n->first << " are in excess!");
  }

  akantu::finalize();

  return EXIT_SUCCESS;
}
