/**
 * @file   mesh_partition_scotch.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  implementation of the MeshPartitionScotch class
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
#include "mesh_partition_scotch.hh"
#include "aka_common.hh"
#include "aka_random_generator.hh"
#include "aka_static_if.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <cstdio>
#include <fstream>
/* -------------------------------------------------------------------------- */

#if !defined(AKANTU_USE_PTSCOTCH)
#ifndef AKANTU_SCOTCH_NO_EXTERN
extern "C" {
#endif // AKANTU_SCOTCH_NO_EXTERN
#include <scotch.h>
#ifndef AKANTU_SCOTCH_NO_EXTERN
}
#endif // AKANTU_SCOTCH_NO_EXTERN
#else  // AKANTU_USE_PTSCOTCH
#include <ptscotch.h>
#endif // AKANTU_USE_PTSCOTCH

namespace akantu {

namespace {
  constexpr int scotch_version = int(SCOTCH_VERSION);
}

/* -------------------------------------------------------------------------- */
MeshPartitionScotch::MeshPartitionScotch(const Mesh & mesh,
                                         UInt spatial_dimension, const ID & id,
                                         const MemoryID & memory_id)
    : MeshPartition(mesh, spatial_dimension, id, memory_id) {
  AKANTU_DEBUG_IN();

  // check if the akantu types and Scotch one are consistent
  static_assert(
      sizeof(Int) == sizeof(SCOTCH_Num),
      "The integer type of Akantu does not match the one from Scotch");

  static_if(aka::bool_constant<scotch_version >= 6>{})
      .then([](auto && y) { SCOTCH_randomSeed(y); })
      .else_([](auto && y) { srandom(y); })(
          std::forward<UInt>(RandomGenerator<UInt>::seed()));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
static SCOTCH_Mesh * createMesh(const Mesh & mesh) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_nodes = mesh.getNbNodes();

  UInt total_nb_element = 0;
  UInt nb_edge = 0;

  for (auto & type : mesh.elementTypes(spatial_dimension)) {
    UInt nb_element = mesh.getNbElement(type);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    total_nb_element += nb_element;
    nb_edge += nb_element * nb_nodes_per_element;
  }

  SCOTCH_Num vnodbas = 0;
  SCOTCH_Num vnodnbr = nb_nodes;

  SCOTCH_Num velmbas = vnodnbr;
  SCOTCH_Num velmnbr = total_nb_element;

  auto * verttab = new SCOTCH_Num[vnodnbr + velmnbr + 1];
  SCOTCH_Num * vendtab = verttab + 1;

  SCOTCH_Num * velotab = nullptr;
  SCOTCH_Num * vnlotab = nullptr;
  SCOTCH_Num * vlbltab = nullptr;

  memset(verttab, 0, (vnodnbr + velmnbr + 1) * sizeof(SCOTCH_Num));

  for (auto & type : mesh.elementTypes(spatial_dimension)) {
    if (Mesh::getSpatialDimension(type) != spatial_dimension)
      continue;

    UInt nb_element = mesh.getNbElement(type);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    const Array<UInt> & connectivity = mesh.getConnectivity(type, _not_ghost);

    /// count number of occurrence of each node
    for (UInt el = 0; el < nb_element; ++el) {
      UInt * conn_val = connectivity.storage() + el * nb_nodes_per_element;
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
        verttab[*(conn_val++)]++;
      }
    }
  }

  /// convert the occurrence array in a csr one
  for (UInt i = 1; i < nb_nodes; ++i)
    verttab[i] += verttab[i - 1];
  for (UInt i = nb_nodes; i > 0; --i)
    verttab[i] = verttab[i - 1];
  verttab[0] = 0;

  /// rearrange element to get the node-element list
  SCOTCH_Num edgenbr = verttab[vnodnbr] + nb_edge;
  auto * edgetab = new SCOTCH_Num[edgenbr];

  UInt linearized_el = 0;

  for (auto & type : mesh.elementTypes(spatial_dimension)) {
    UInt nb_element = mesh.getNbElement(type);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    const Array<UInt> & connectivity = mesh.getConnectivity(type, _not_ghost);

    for (UInt el = 0; el < nb_element; ++el, ++linearized_el) {
      UInt * conn_val = connectivity.storage() + el * nb_nodes_per_element;
      for (UInt n = 0; n < nb_nodes_per_element; ++n)
        edgetab[verttab[*(conn_val++)]++] = linearized_el + velmbas;
    }
  }

  for (UInt i = nb_nodes; i > 0; --i)
    verttab[i] = verttab[i - 1];
  verttab[0] = 0;

  SCOTCH_Num * verttab_tmp = verttab + vnodnbr + 1;
  SCOTCH_Num * edgetab_tmp = edgetab + verttab[vnodnbr];

  for (auto & type : mesh.elementTypes(spatial_dimension)) {
    UInt nb_element = mesh.getNbElement(type);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    const Array<UInt> & connectivity = mesh.getConnectivity(type, _not_ghost);

    for (UInt el = 0; el < nb_element; ++el) {
      *verttab_tmp = *(verttab_tmp - 1) + nb_nodes_per_element;
      verttab_tmp++;
      UInt * conn = connectivity.storage() + el * nb_nodes_per_element;
      for (UInt i = 0; i < nb_nodes_per_element; ++i) {
        *(edgetab_tmp++) = *(conn++) + vnodbas;
      }
    }
  }

  auto * meshptr = new SCOTCH_Mesh;

  SCOTCH_meshInit(meshptr);

  SCOTCH_meshBuild(meshptr, velmbas, vnodbas, velmnbr, vnodnbr, verttab,
                   vendtab, velotab, vnlotab, vlbltab, edgenbr, edgetab);

  /// Check the mesh
  AKANTU_DEBUG_ASSERT(SCOTCH_meshCheck(meshptr) == 0,
                      "Scotch mesh is not consistent");

#ifndef AKANTU_NDEBUG
  if (AKANTU_DEBUG_TEST(dblDump)) {
    /// save initial graph
    FILE * fmesh = fopen("ScotchMesh.msh", "w");
    SCOTCH_meshSave(meshptr, fmesh);
    fclose(fmesh);

    /// write geometry file
    std::ofstream fgeominit;
    fgeominit.open("ScotchMesh.xyz");
    fgeominit << spatial_dimension << std::endl << nb_nodes << std::endl;

    const Array<Real> & nodes = mesh.getNodes();
    Real * nodes_val = nodes.storage();
    for (UInt i = 0; i < nb_nodes; ++i) {
      fgeominit << i << " ";
      for (UInt s = 0; s < spatial_dimension; ++s)
        fgeominit << *(nodes_val++) << " ";
      fgeominit << std::endl;
      ;
    }
    fgeominit.close();
  }
#endif

  AKANTU_DEBUG_OUT();
  return meshptr;
}

/* -------------------------------------------------------------------------- */
static void destroyMesh(SCOTCH_Mesh * meshptr) {
  AKANTU_DEBUG_IN();

  SCOTCH_Num velmbas, vnodbas, vnodnbr, velmnbr, *verttab, *vendtab, *velotab,
      *vnlotab, *vlbltab, edgenbr, *edgetab, degrptr;

  SCOTCH_meshData(meshptr, &velmbas, &vnodbas, &velmnbr, &vnodnbr, &verttab,
                  &vendtab, &velotab, &vnlotab, &vlbltab, &edgenbr, &edgetab,
                  &degrptr);

  delete[] verttab;
  delete[] edgetab;

  SCOTCH_meshExit(meshptr);

  delete meshptr;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshPartitionScotch::partitionate(
    UInt nb_part,
    std::function<Int(const Element &, const Element &)> edge_load_func,
    std::function<Int(const Element &)> vertex_load_func) {
  AKANTU_DEBUG_IN();

  nb_partitions = nb_part;

  if (mesh.isPeriodic()) {
    tweakConnectivity();
  }

  AKANTU_DEBUG_INFO("Partitioning the mesh " << mesh.getID() << " in "
                                             << nb_part << " parts.");

  Array<Int> dxadj;
  Array<Int> dadjncy;
  Array<Int> edge_loads;
  Array<Int> vertex_loads;
  buildDualGraph(dxadj, dadjncy, edge_loads, edge_load_func, vertex_loads,
                 vertex_load_func);

  /// variables that will hold our structures in scotch format
  SCOTCH_Graph scotch_graph;
  SCOTCH_Strat scotch_strat;

  /// description number and arrays for struct mesh for scotch
  SCOTCH_Num baseval = 0; // base numbering for element and
  // nodes (0 -> C , 1 -> fortran)
  SCOTCH_Num vertnbr = dxadj.size() - 1; // number of vertexes
  SCOTCH_Num * parttab;                  // array of partitions
  SCOTCH_Num edgenbr = dxadj(vertnbr);   // twice  the number  of "edges"
  //(an "edge" bounds two nodes)
  SCOTCH_Num * verttab = dxadj.storage(); // array of start indices in edgetab
  SCOTCH_Num * vendtab = nullptr; // array of after-last indices in edgetab
  SCOTCH_Num * velotab =
      vertex_loads.storage(); // integer  load  associated with
  // every vertex ( optional )
  SCOTCH_Num * edlotab = edge_loads.storage(); // integer  load  associated with
  // every edge ( optional )
  SCOTCH_Num * edgetab = dadjncy.storage(); // adjacency array of every vertex
  SCOTCH_Num * vlbltab = nullptr;           // vertex label array (optional)

  /// Allocate space for Scotch arrays
  parttab = new SCOTCH_Num[vertnbr];

  /// Initialize the strategy structure
  SCOTCH_stratInit(&scotch_strat);

  /// Initialize the graph structure
  SCOTCH_graphInit(&scotch_graph);

  /// Build the graph from the adjacency arrays
  SCOTCH_graphBuild(&scotch_graph, baseval, vertnbr, verttab, vendtab, velotab,
                    vlbltab, edgenbr, edgetab, edlotab);

#ifndef AKANTU_NDEBUG
  if (AKANTU_DEBUG_TEST(dblDump)) {
    /// save initial graph
    FILE * fgraphinit = fopen("GraphIniFile.grf", "w");
    SCOTCH_graphSave(&scotch_graph, fgraphinit);
    fclose(fgraphinit);

    /// write geometry file
    std::ofstream fgeominit;
    fgeominit.open("GeomIniFile.xyz");
    fgeominit << spatial_dimension << std::endl << vertnbr << std::endl;

    const Array<Real> & nodes = mesh.getNodes();

    auto nodes_it = nodes.begin(spatial_dimension);

    UInt out_linerized_el = 0;
    for (auto & type :
         mesh.elementTypes(spatial_dimension, _not_ghost, _ek_not_defined)) {
      UInt nb_element = mesh.getNbElement(type);
      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
      const Array<UInt> & connectivity = mesh.getConnectivity(type);

      Vector<Real> mid(spatial_dimension);
      for (UInt el = 0; el < nb_element; ++el) {
        mid.set(0.);
        for (UInt n = 0; n < nb_nodes_per_element; ++n) {
          UInt node = connectivity.storage()[nb_nodes_per_element * el + n];
          mid += Vector<Real>(nodes_it[node]);
        }
        mid /= nb_nodes_per_element;

        fgeominit << out_linerized_el++ << " ";
        for (UInt s = 0; s < spatial_dimension; ++s)
          fgeominit << mid[s] << " ";

        fgeominit << std::endl;
        ;
      }
    }
    fgeominit.close();
  }
#endif

  /// Check the graph
  AKANTU_DEBUG_ASSERT(SCOTCH_graphCheck(&scotch_graph) == 0,
                      "Graph to partition is not consistent");

  /// Partition the mesh
  SCOTCH_graphPart(&scotch_graph, nb_part, &scotch_strat, parttab);

  /// Check the graph
  AKANTU_DEBUG_ASSERT(SCOTCH_graphCheck(&scotch_graph) == 0,
                      "Partitioned graph is not consistent");

#ifndef AKANTU_NDEBUG
  if (AKANTU_DEBUG_TEST(dblDump)) {
    /// save the partitioned graph
    FILE * fgraph = fopen("GraphFile.grf", "w");
    SCOTCH_graphSave(&scotch_graph, fgraph);
    fclose(fgraph);

    /// save the partition map
    std::ofstream fmap;
    fmap.open("MapFile.map");
    fmap << vertnbr << std::endl;
    for (Int i = 0; i < vertnbr; i++)
      fmap << i << "    " << parttab[i] << std::endl;
    fmap.close();
  }
#endif

  /// free the scotch data structures
  SCOTCH_stratExit(&scotch_strat);
  SCOTCH_graphFree(&scotch_graph);
  SCOTCH_graphExit(&scotch_graph);

  fillPartitionInformation(mesh, parttab);

  delete[] parttab;

  if (mesh.isPeriodic()) {
    restoreConnectivity();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshPartitionScotch::reorder() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Reordering the mesh " << mesh.getID());
  SCOTCH_Mesh * scotch_mesh = createMesh(mesh);

  UInt nb_nodes = mesh.getNbNodes();

  SCOTCH_Strat scotch_strat;
  // SCOTCH_Ordering scotch_order;

  auto * permtab = new SCOTCH_Num[nb_nodes];
  SCOTCH_Num * peritab = nullptr;
  SCOTCH_Num cblknbr = 0;
  SCOTCH_Num * rangtab = nullptr;
  SCOTCH_Num * treetab = nullptr;

  /// Initialize the strategy structure
  SCOTCH_stratInit(&scotch_strat);

  SCOTCH_Graph scotch_graph;

  SCOTCH_graphInit(&scotch_graph);
  SCOTCH_meshGraph(scotch_mesh, &scotch_graph);

#ifndef AKANTU_NDEBUG
  if (AKANTU_DEBUG_TEST(dblDump)) {
    FILE * fgraphinit = fopen("ScotchMesh.grf", "w");
    SCOTCH_graphSave(&scotch_graph, fgraphinit);
    fclose(fgraphinit);
  }
#endif

  /// Check the graph
  // AKANTU_DEBUG_ASSERT(SCOTCH_graphCheck(&scotch_graph) == 0,
  //		      "Mesh to Graph is not consistent");

  SCOTCH_graphOrder(&scotch_graph, &scotch_strat, permtab, peritab, &cblknbr,
                    rangtab, treetab);

  SCOTCH_graphExit(&scotch_graph);
  SCOTCH_stratExit(&scotch_strat);
  destroyMesh(scotch_mesh);

  /// Renumbering
  UInt spatial_dimension = mesh.getSpatialDimension();

  for (auto gt : ghost_types) {
    for (auto & type : mesh.elementTypes(_ghost_type = gt)) {
      UInt nb_element = mesh.getNbElement(type, gt);
      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

      const Array<UInt> & connectivity = mesh.getConnectivity(type, gt);

      UInt * conn = connectivity.storage();
      for (UInt el = 0; el < nb_element * nb_nodes_per_element; ++el, ++conn) {
        *conn = permtab[*conn];
      }
    }
  }

  /// \todo think of a in-place way to do it
  auto * new_coordinates = new Real[spatial_dimension * nb_nodes];
  Real * old_coordinates = mesh.getNodes().storage();
  for (UInt i = 0; i < nb_nodes; ++i) {
    memcpy(new_coordinates + permtab[i] * spatial_dimension,
           old_coordinates + i * spatial_dimension,
           spatial_dimension * sizeof(Real));
  }
  memcpy(old_coordinates, new_coordinates,
         nb_nodes * spatial_dimension * sizeof(Real));
  delete[] new_coordinates;

  delete[] permtab;

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
