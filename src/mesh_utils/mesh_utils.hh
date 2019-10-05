/**
 * @file   mesh_utils.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Dec 03 2017
 *
 * @brief  All mesh utils necessary for various tasks
 *
 * @section LICENSE
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
#include "aka_csr.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MESH_UTILS_HH__
#define __AKANTU_MESH_UTILS_HH__

namespace akantu {

class MeshUtils {
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// build a CSR<Element> that contains for each node the list of connected
  /// elements of a given spatial dimension
  static void buildNode2Elements(const Mesh & mesh, CSR<Element> & node_to_elem,
                                 UInt spatial_dimension = _all_dimensions);

  /// build a CSR<UInt> that contains for each node the number of
  /// the connected elements of a given ElementType
  static void
  buildNode2ElementsElementTypeMap(const Mesh & mesh, CSR<UInt> & node_to_elem,
                                   const ElementType & type,
                                   const GhostType & ghost_type = _not_ghost);

  /// build the facets elements on the boundaries of a mesh
  static void buildFacets(Mesh & mesh);

  /// build all the facets elements: boundary and internals and store them in
  /// the mesh_facets for element of dimension from_dimension to to_dimension
  static void buildAllFacets(const Mesh & mesh, Mesh & mesh_facets,
                             UInt from_dimension, UInt to_dimension);

  /// build all the facets elements: boundary and internals and store them in
  /// the mesh_facets
  static void buildAllFacets(const Mesh & mesh, Mesh & mesh_facets,
                             UInt to_dimension = 0);

  /// build facets for a given spatial dimension
  static void buildFacetsDimension(const Mesh & mesh, Mesh & mesh_facets,
                                   bool boundary_only, UInt dimension);

  /// take the local_connectivity array as the array of local and ghost
  /// connectivity, renumber the nodes and set the connectivity of the mesh
  static void renumberMeshNodes(Mesh & mesh, Array<UInt> & local_connectivities,
                                UInt nb_local_element, UInt nb_ghost_element,
                                ElementType type, Array<UInt> & old_nodes);

  /// compute pbc pair for a given direction
  static void computePBCMap(const Mesh & mymesh, const UInt dir,
                            std::map<UInt, UInt> & pbc_pair);
  /// compute pbc pair for a surface pair
  static void computePBCMap(const Mesh & mymesh,
                            const std::pair<ID, ID> & surface_pair,
                            std::map<UInt, UInt> & pbc_pair);

  /// remove not connected nodes /!\ this functions renumbers the nodes.
  static void purifyMesh(Mesh & mesh);

#if defined(AKANTU_COHESIVE_ELEMENT)
  /// function to insert cohesive elements on the selected facets
  /// @return number of facets that have been doubled
  static UInt
  insertCohesiveElements(Mesh & mesh, Mesh & mesh_facets,
                         const ElementTypeMapArray<bool> & facet_insertion,
                         Array<UInt> & doubled_nodes,
                         Array<Element> & new_elements,
                         bool only_double_facets);
#endif

  /// fill the subelement to element and the elements to subelements data
  static void fillElementToSubElementsData(Mesh & mesh);

  /// flip facets based on global connectivity
  static void flipFacets(Mesh & mesh_facets,
                         const ElementTypeMapArray<UInt> & global_connectivity,
                         GhostType gt_facet);

  /// provide list of elements around a node and check if a given
  /// facet is reached
  template <bool third_dim_points>
  static bool findElementsAroundSubfacet(
      const Mesh & mesh_facets, const Element & starting_element,
      const Element & end_facet, const Vector<UInt> & subfacet_connectivity,
      std::vector<Element> & elem_list, std::vector<Element> & facet_list,
      std::vector<Element> * subfacet_list = nullptr);

  /// function to check if a node belongs to a given element
  static inline bool hasElement(const Vector<UInt> & nodes_element,
                                const Vector<UInt> & nodes);

  /// reset facet_to_double arrays in the Mesh
  static void resetFacetToDouble(Mesh & mesh_facets);

private:
  /// match pairs that are on the associated pbc's
  static void matchPBCPairs(const Mesh & mymesh, const UInt dir,
                            Array<UInt> & selected_left,
                            Array<UInt> & selected_right,
                            std::map<UInt, UInt> & pbc_pair);

  /// function used by all the renumbering functions
  static void
  renumberNodesInConnectivity(Array<UInt> & list_nodes, UInt nb_nodes,
                              std::map<UInt, UInt> & renumbering_map);

  /// update facet_to_subfacet
  static void updateFacetToSubfacet(Mesh & mesh_facets,
                                    ElementType type_subfacet,
                                    GhostType gt_subfacet, bool facet_mode);

  /// update subfacet_to_facet
  static void updateSubfacetToFacet(Mesh & mesh_facets,
                                    ElementType type_subfacet,
                                    GhostType gt_subfacet, bool facet_mode);

  /// function to double a given facet and update the list of doubled
  /// nodes
  static void doubleFacet(Mesh & mesh, Mesh & mesh_facets, UInt facet_dimension,
                          Array<UInt> & doubled_nodes, bool facet_mode);

  /// function to double a subfacet given start and end index for
  /// local facet_to_subfacet vector, and update the list of doubled
  /// nodes
  template <UInt spatial_dimension>
  static void doubleSubfacet(Mesh & mesh, Mesh & mesh_facets,
                             Array<UInt> & doubled_nodes);

  /// double a node
  static void doubleNodes(Mesh & mesh, const std::vector<UInt> & old_nodes,
                          Array<UInt> & doubled_nodes);

  /// fill facet_to_double array in the mesh
  /// returns the number of facets to be doubled
  static UInt
  updateFacetToDouble(Mesh & mesh_facets,
                      const ElementTypeMapArray<bool> & facet_insertion);

  /// find subfacets to be doubled
  template <bool subsubfacet_mode>
  static void findSubfacetToDouble(Mesh & mesh_facets);

  /// double facets (points) in 1D
  static void doublePointFacet(Mesh & mesh, Mesh & mesh_facets,
                               Array<UInt> & doubled_nodes);

#if defined(AKANTU_COHESIVE_ELEMENT)
  /// update cohesive element data
  static void updateCohesiveData(Mesh & mesh, Mesh & mesh_facets,
                                 Array<Element> & new_elements);
#endif

  /// update elemental connectivity after doubling a node
  inline static void updateElementalConnectivity(
      Mesh & mesh, UInt old_node, UInt new_node,
      const std::vector<Element> & element_list,
      const std::vector<Element> * facet_list = nullptr);

  /// double middle nodes if facets are _segment_3
  template <bool third_dim_segments>
  static void updateQuadraticSegments(Mesh & mesh, Mesh & mesh_facets,
                                      ElementType type_facet,
                                      GhostType gt_facet,
                                      Array<UInt> & doubled_nodes);

  /// remove elements on a vector
  inline static bool
  removeElementsInVector(const std::vector<Element> & elem_to_remove,
                         std::vector<Element> & elem_list);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "mesh_utils_inline_impl.cc"

#endif /* __AKANTU_MESH_UTILS_HH__ */
