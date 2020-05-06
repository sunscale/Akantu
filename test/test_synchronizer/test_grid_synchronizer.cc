/**
 * @file   test_grid_synchronizer.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 01 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  test the GridSynchronizer object
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
#include "aka_common.hh"
#include "aka_grid_dynamic.hh"
#include "grid_synchronizer.hh"
#include "mesh.hh"
#include "mesh_partition.hh"
#include "synchronizer_registry.hh"
#include "test_data_accessor.hh"
#ifdef AKANTU_USE_IOHELPER
#include "io_helper.hh"
#endif // AKANTU_USE_IOHELPER

using namespace akantu;

const UInt spatial_dimension = 2;

typedef std::map<std::pair<Element, Element>, Real> pair_list;

#include "test_grid_tools.hh"

static void
updatePairList(const ElementTypeMapArray<Real> & barycenter,
               const SpatialGrid<Element> & grid, Real radius,
               pair_list & neighbors,
               neighbors_map_t<spatial_dimension>::type & neighbors_map) {
  AKANTU_DEBUG_IN();

  GhostType ghost_type = _not_ghost;

  Element e;
  e.ghost_type = ghost_type;

  // generate the pair of neighbor depending of the cell_list
  ElementTypeMapArray<Real>::type_iterator it =
      barycenter.firstType(_all_dimensions, ghost_type);
  ElementTypeMapArray<Real>::type_iterator last_type =
      barycenter.lastType(0, ghost_type);
  for (; it != last_type; ++it) {
    // loop over quad points

    e.type = *it;
    e.element = 0;

    const Array<Real> & barycenter_vect = barycenter(*it, ghost_type);
    UInt sp = barycenter_vect.getNbComponent();

    Array<Real>::const_iterator<Vector<Real>> bary = barycenter_vect.begin(sp);
    Array<Real>::const_iterator<Vector<Real>> bary_end =
        barycenter_vect.end(sp);

    for (; bary != bary_end; ++bary, e.element++) {
#if !defined(AKANTU_NDEBUG)
      Point<spatial_dimension> pt1(*bary);
#endif

      SpatialGrid<Element>::CellID cell_id = grid.getCellID(*bary);
      SpatialGrid<Element>::neighbor_cells_iterator first_neigh_cell =
          grid.beginNeighborCells(cell_id);
      SpatialGrid<Element>::neighbor_cells_iterator last_neigh_cell =
          grid.endNeighborCells(cell_id);

      // loop over neighbors cells of the one containing the current element
      for (; first_neigh_cell != last_neigh_cell; ++first_neigh_cell) {
        SpatialGrid<Element>::Cell::const_iterator first_neigh_el =
            grid.beginCell(*first_neigh_cell);
        SpatialGrid<Element>::Cell::const_iterator last_neigh_el =
            grid.endCell(*first_neigh_cell);

        // loop over the quadrature point in the current cell of the cell list
        for (; first_neigh_el != last_neigh_el; ++first_neigh_el) {
          const Element & elem = *first_neigh_el;

          Array<Real>::const_iterator<Vector<Real>> neigh_it =
              barycenter(elem.type, elem.ghost_type).begin(sp);

          const Vector<Real> & neigh_bary = neigh_it[elem.element];

          Real distance = bary->distance(neigh_bary);
          if (distance <= radius) {
#if !defined(AKANTU_NDEBUG)
            Point<spatial_dimension> pt2(neigh_bary);
            neighbors_map[pt1].push_back(pt2);
#endif
            std::pair<Element, Element> pair = std::make_pair(e, elem);
            pair_list::iterator p = neighbors.find(pair);
            if (p != neighbors.end()) {
              AKANTU_ERROR("Pair already registered ["
                           << e << " " << elem << "] -> " << p->second << " "
                           << distance);
            } else {
              neighbors[pair] = distance;
            }
          }
        }
      }
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  akantu::initialize(argc, argv);

  Real radius = 0.001;

  Mesh mesh(spatial_dimension);

  const auto & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();
  ElementSynchronizer * dist = NULL;

  if (prank == 0) {
    mesh.read("bar.msh");
    MeshPartition * partition =
        new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
    dist =
        ElementSynchronizer::createDistributedSynchronizerMesh(mesh, partition);
    delete partition;
  } else {
    dist = ElementSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }

  mesh.computeBoundingBox();

  const Vector<Real> & lower_bounds = mesh.getLowerBounds();
  const Vector<Real> & upper_bounds = mesh.getUpperBounds();

  Vector<Real> center = 0.5 * (upper_bounds + lower_bounds);

  Vector<Real> spacing(spatial_dimension);
  for (UInt i = 0; i < spatial_dimension; ++i) {
    spacing[i] = radius * 1.2;
  }

  SpatialGrid<Element> grid(spatial_dimension, spacing, center);

  GhostType ghost_type = _not_ghost;

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);

  ElementTypeMapArray<Real> barycenters("", "");
  mesh.initElementTypeMapArray(barycenters, spatial_dimension,
                               spatial_dimension);

  Element e;
  e.ghost_type = ghost_type;

  for (; it != last_type; ++it) {
    UInt nb_element = mesh.getNbElement(*it, ghost_type);
    e.type = *it;
    Array<Real> & barycenter = barycenters(*it, ghost_type);
    barycenter.resize(nb_element);

    Array<Real>::iterator<Vector<Real>> bary_it =
        barycenter.begin(spatial_dimension);
    for (UInt elem = 0; elem < nb_element; ++elem) {
      mesh.getBarycenter(elem, *it, bary_it->storage(), ghost_type);
      e.element = elem;
      grid.insert(e, *bary_it);
      ++bary_it;
    }
  }

  std::stringstream sstr;
  sstr << "mesh_" << prank << ".msh";
  mesh.write(sstr.str());

  Mesh grid_mesh(spatial_dimension, "grid_mesh", 0);
  std::stringstream sstr_grid;
  sstr_grid << "grid_mesh_" << prank << ".msh";
  grid.saveAsMesh(grid_mesh);
  grid_mesh.write(sstr_grid.str());

  std::cout << "Pouet 1" << std::endl;

  AKANTU_DEBUG_INFO("Creating TestAccessor");
  TestAccessor test_accessor(mesh, barycenters);
  SynchronizerRegistry synch_registry(test_accessor);

  GridSynchronizer * grid_communicator =
      GridSynchronizer::createGridSynchronizer(mesh, grid);

  std::cout << "Pouet 2" << std::endl;

  ghost_type = _ghost;

  it = mesh.firstType(spatial_dimension, ghost_type);
  last_type = mesh.lastType(spatial_dimension, ghost_type);
  e.ghost_type = ghost_type;
  for (; it != last_type; ++it) {
    UInt nb_element = mesh.getNbElement(*it, ghost_type);
    e.type = *it;
    Array<Real> & barycenter = barycenters(*it, ghost_type);
    barycenter.resize(nb_element);

    Array<Real>::iterator<Vector<Real>> bary_it =
        barycenter.begin(spatial_dimension);
    for (UInt elem = 0; elem < nb_element; ++elem) {
      mesh.getBarycenter(elem, *it, bary_it->storage(), ghost_type);
      e.element = elem;
      grid.insert(e, *bary_it);
      ++bary_it;
    }
  }

  Mesh grid_mesh_ghost(spatial_dimension, "grid_mesh_ghost", 0);
  std::stringstream sstr_gridg;
  sstr_gridg << "grid_mesh_ghost_" << prank << ".msh";
  grid.saveAsMesh(grid_mesh_ghost);
  grid_mesh_ghost.write(sstr_gridg.str());

  std::cout << "Pouet 3" << std::endl;

  neighbors_map_t<spatial_dimension>::type neighbors_map;
  pair_list neighbors;

  updatePairList(barycenters, grid, radius, neighbors, neighbors_map);
  pair_list::iterator nit = neighbors.begin();
  pair_list::iterator nend = neighbors.end();

  std::stringstream sstrp;
  sstrp << "pairs_" << prank;
  std::ofstream fout(sstrp.str().c_str());
  for (; nit != nend; ++nit) {
    fout << "[" << nit->first.first << "," << nit->first.second << "] -> "
         << nit->second << std::endl;
  }

  std::string file = "neighbors_ref";
  std::stringstream sstrf;
  sstrf << file << "_" << psize << "_" << prank;

  file = sstrf.str();

  std::ofstream nout;
  nout.open(file.c_str());
  neighbors_map_t<spatial_dimension>::type::iterator it_n =
      neighbors_map.begin();
  neighbors_map_t<spatial_dimension>::type::iterator end_n =
      neighbors_map.end();
  for (; it_n != end_n; ++it_n) {
    std::sort(it_n->second.begin(), it_n->second.end());

    std::vector<Point<spatial_dimension>>::iterator it_v = it_n->second.begin();
    std::vector<Point<spatial_dimension>>::iterator end_v = it_n->second.end();

    nout << "####" << std::endl;
    nout << it_n->second.size() << std::endl;
    nout << it_n->first << std::endl;
    nout << "#" << std::endl;
    for (; it_v != end_v; ++it_v) {
      nout << *it_v << std::endl;
    }
  }

  fout.close();

  synch_registry.registerSynchronizer(*dist, SynchronizationTag::_smm_mass);

  synch_registry.registerSynchronizer(*grid_communicator,
                                      SynchronizationTag::_test);

  AKANTU_DEBUG_INFO("Synchronizing tag on Dist");
  synch_registry.synchronize(SynchronizationTag::_smm_mass);

  AKANTU_DEBUG_INFO("Synchronizing tag on Grid");
  synch_registry.synchronize(SynchronizationTag::_test);

  delete grid_communicator;
  delete dist;
  akantu::finalize();

  return EXIT_SUCCESS;
}
