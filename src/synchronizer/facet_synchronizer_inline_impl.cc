/**
 * @file   facet_synchronizer_inline_impl.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 *
 * @brief  facet synchronizer inline implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
// template<GhostType ghost_facets>
// inline void FacetSynchronizer::getFacetGlobalConnectivity(const DistributedSynchronizer & distributed_synchronizer,
//                                                        const ElementTypeMapArray<UInt> & rank_to_facet,
//                                                        const Array<Element> * elements,
//                                                        Array<ElementTypeMapArray<UInt> *> & connectivity,
//                                                        Array<ElementTypeMapArray<UInt> *> & facets) {
//   AKANTU_DEBUG_IN();

//   UInt spatial_dimension = mesh.getSpatialDimension();

//   /// init facet check tracking
//   ElementTypeMapArray<bool> facet_checked("facet_checked", id);

//   mesh.initElementTypeMapArray(facet_checked, 1, spatial_dimension - 1);

//   Mesh::type_iterator first = mesh.firstType(spatial_dimension - 1, ghost_facets);
//   Mesh::type_iterator last  = mesh.lastType(spatial_dimension - 1, ghost_facets);

//   for (; first != last; ++first) {
//     ElementType type = *first;
//     Array<bool> & f_checked = facet_checked(type, ghost_facets);
//     UInt nb_facet = mesh.getNbElement(type, ghost_facets);
//     f_checked.resize(nb_facet);
//   }

//   const Array<UInt> & nodes_global_ids =
//     distributed_synchronizer.mesh.getGlobalNodesIds();

//   /// loop on every processor
//   for (UInt p = 0; p < nb_proc; ++p) {
//     if (p == rank) continue;

//     /// reset facet check
//     Mesh::type_iterator first = mesh.firstType(spatial_dimension - 1, ghost_facets);
//     Mesh::type_iterator last  = mesh.lastType(spatial_dimension - 1, ghost_facets);

//     for (; first != last; ++first) {
//       ElementType type = *first;
//       Array<bool> & f_checked = facet_checked(type, ghost_facets);
//       f_checked.clear();
//     }

//     ElementTypeMapArray<UInt> & global_conn = (*connectivity(p));
//     const Array<Element> & elem = elements[p];
//     ElementTypeMapArray<UInt> & facet_list = (*facets(p));

//     UInt nb_element = elem.getSize();

//     /// loop on every send/recv element
//     for (UInt el = 0; el < nb_element; ++el) {
//       ElementType type = elem(el).type;
//       GhostType gt = elem(el).ghost_type;
//       UInt el_index = elem(el).element;

//       const Array<Element> & facet_to_element =
//      mesh.getSubelementToElement(type, gt);
//       UInt nb_facets_per_element = Mesh::getNbFacetsPerElement(type);
//       ElementType facet_type = Mesh::getFacetType(type);
//       UInt nb_nodes_per_facet = Mesh::getNbNodesPerElement(facet_type);
//       Vector<UInt> conn_tmp(nb_nodes_per_facet);

//       /// loop on every facet of the element
//       for (UInt f = 0; f < nb_facets_per_element; ++f) {

//      const Element & facet = facet_to_element(el_index, f);
//      if (facet == ElementNull) continue;
//      UInt facet_index = facet.element;
//      GhostType facet_gt = facet.ghost_type;

//      const Array<UInt> & t_to_f = rank_to_facet(facet_type, facet_gt);

//      /// exclude not ghost facets, facets assigned to other
//      /// processors
//      if (facet_gt != ghost_facets) continue;
//      if ((facet_gt == _ghost) && (t_to_f(facet_index) != p)) continue;

//      /// exclude facets that have already been added
//      Array<bool> & f_checked = facet_checked(facet_type, facet_gt);
//      if (f_checked(facet_index)) continue;
//      else f_checked(facet_index) = true;

//      /// add facet index
//      Array<UInt> & f_list = facet_list(facet_type, facet_gt);
//      f_list.push_back(facet_index);

//      /// add sorted facet global connectivity
//      const Array<UInt> & conn = mesh.getConnectivity(facet_type, facet_gt);
//      Array<UInt> & g_conn = global_conn(facet_type, facet_gt);

//      for (UInt n = 0; n < nb_nodes_per_facet; ++n)
//        conn_tmp(n) = nodes_global_ids(conn(facet_index, n));

//      std::sort(conn_tmp.storage(), conn_tmp.storage() + nb_nodes_per_facet);

//      g_conn.push_back(conn_tmp);
//       }
//     }
//   }

//   AKANTU_DEBUG_OUT();
// }
