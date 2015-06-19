/**
 * @file   element_class_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Specialization for interface-enriched finite elements
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

#include "element_class.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_igfem_segment_3>::spatial_dimension          = 1;
template<> UInt GeometricalElement<_gt_igfem_segment_3>::nb_nodes_per_element       = 3;
template<> UInt GeometricalElement<_gt_igfem_segment_3>::nb_facets[]                = { 2 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_igfem_triangle_4>::spatial_dimension          = 2;
template<> UInt GeometricalElement<_gt_igfem_triangle_4>::nb_nodes_per_element       = 4;
template<> UInt GeometricalElement<_gt_igfem_triangle_4>::nb_facets[]                = { 3 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_igfem_triangle_5>::spatial_dimension          = 2;
template<> UInt GeometricalElement<_gt_igfem_triangle_5>::nb_nodes_per_element       = 5;
template<> UInt GeometricalElement<_gt_igfem_triangle_5>::nb_facets[]                = { 3 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_igfem_segment_3>::nb_nodes_per_facet[]      = { 1 };
template<> UInt GeometricalElement<_gt_igfem_triangle_4>::nb_nodes_per_facet[]     = { 2, 3 };
template<> UInt GeometricalElement<_gt_igfem_triangle_5>::nb_nodes_per_facet[]     = { 2, 3 };
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_igfem_segment_3>::nb_facet_types      = 1;
template<> UInt GeometricalElement<_gt_igfem_triangle_4>::nb_facet_types     = 2;
template<> UInt GeometricalElement<_gt_igfem_triangle_5>::nb_facet_types     = 2;
/* !!! stored as a matrix nb_facets X nb_nodes_per_facet in COL MAJOR */
/* -------------------------------------------------------------------------- */
/*                                                                                  f1|f2|f3|f4|f5|f6 */
template<> UInt GeometricalElement<_gt_igfem_segment_3>::facet_connectivity_vect[]      = {0, 1};
template<> UInt GeometricalElement<_gt_igfem_triangle_4>::facet_connectivity_vect[]     = {// first type
                                                                                     0, 2,
										     1, 0,
										     // second type
										     1,  
										     2,
										     3};
template<> UInt GeometricalElement<_gt_igfem_triangle_5>::facet_connectivity_vect[]     = {// first type
                                                                                     1,
										     2,
										     // second type
										     0, 2,
										     1, 0,
										     3, 4};

template<> UInt * GeometricalElement<_gt_igfem_segment_3>::facet_connectivity[]      = { &facet_connectivity_vect[0] };
template<> UInt * GeometricalElement<_gt_igfem_triangle_4>::facet_connectivity[]     = { &facet_connectivity_vect[0], &facet_connectivity_vect[4] };
template<> UInt * GeometricalElement<_gt_igfem_triangle_5>::facet_connectivity[]     = { &facet_connectivity_vect[0], &facet_connectivity_vect[2] };;

/* -------------------------------------------------------------------------- */

__END_AKANTU__
