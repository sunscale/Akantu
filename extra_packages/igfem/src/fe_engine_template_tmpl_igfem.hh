/**
 * @file   shape_igfem_inline_impl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @brief  ShapeIGFEM inline implementation
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "integrator_gauss_igfem.hh"
#include "shape_igfem.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_FE_ENGINE_TEMPLATE_TMPL_IGFEM_HH_
#define AKANTU_FE_ENGINE_TEMPLATE_TMPL_IGFEM_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
/* compatibility functions */
/* -------------------------------------------------------------------------- */

template <>
inline void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_igfem,
                             DefaultIntegrationOrderFunctor>::
    initShapeFunctions(const Array<Real> & nodes,
                       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it =
      mesh.firstType(element_dimension, ghost_type, _ek_igfem);
  Mesh::type_iterator end =
      mesh.lastType(element_dimension, ghost_type, _ek_igfem);
  for (; it != end; ++it) {
    ElementType type = *it;
    integrator.initIntegrator(nodes, type, ghost_type);

#define INIT(_type)                                                            \
  do {                                                                         \
    const Matrix<Real> & all_quads =                                           \
        integrator.getIntegrationPoints<_type>(ghost_type);                    \
    const Matrix<Real> & quads_1 = integrator.getIntegrationPoints<            \
        ElementClassProperty<_type>::sub_element_type_1>(ghost_type);          \
    const Matrix<Real> & quads_2 = integrator.getIntegrationPoints<            \
        ElementClassProperty<_type>::sub_element_type_2>(ghost_type);          \
    shape_functions.initShapeFunctions(nodes, all_quads, quads_1, quads_2,     \
                                       _type, ghost_type);                     \
  } while (0)

    AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(INIT);
#undef INIT
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <>
inline void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_igfem,
                             DefaultIntegrationOrderFunctor>::
    computeIntegrationPointsCoordinates(
        Array<Real> & quadrature_points_coordinates, ElementType type,
        GhostType ghost_type,
        const Array<UInt> & filter_elements) const {

  const Array<Real> & nodes_coordinates = mesh.getNodes();
  UInt spatial_dimension = mesh.getSpatialDimension();
  /// create an array with the nodal coordinates that need to be
  /// interpolated. The nodal coordinates of the enriched nodes need
  /// to be set to zero, because they represent the enrichment of the
  /// position field, and the enrichments for this field are all zero!
  /// There is no gap in the mesh!
  Array<Real> igfem_nodes(nodes_coordinates.getSize(), spatial_dimension);
  shape_functions.extractValuesAtStandardNodes(nodes_coordinates, igfem_nodes,
                                               ghost_type);

  interpolateOnIntegrationPoints(igfem_nodes, quadrature_points_coordinates,
                                 spatial_dimension, type, ghost_type,
                                 filter_elements);
}

/* -------------------------------------------------------------------------- */
template <>
inline void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_igfem,
                             DefaultIntegrationOrderFunctor>::
    computeIntegrationPointsCoordinates(
        ElementTypeMapArray<Real> & quadrature_points_coordinates,
        const ElementTypeMapArray<UInt> * filter_elements) const {

  const Array<Real> & nodes_coordinates = mesh.getNodes();
  UInt spatial_dimension = mesh.getSpatialDimension();
  /// create an array with the nodal coordinates that need to be
  /// interpolated. The nodal coordinates of the enriched nodes need
  /// to be set to zero, because they represent the enrichment of the
  /// position field, and the enrichments for this field are all zero!
  /// There is no gap in the mesh!
  Array<Real> igfem_nodes(nodes_coordinates.getSize(), spatial_dimension);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;
    shape_functions.extractValuesAtStandardNodes(nodes_coordinates, igfem_nodes,
                                                 ghost_type);
  }

  interpolateOnIntegrationPoints(igfem_nodes, quadrature_points_coordinates,
                                 filter_elements);
}

} // namespace akantu

#endif /* AKANTU_FE_ENGINE_TEMPLATE_TMPL_IGFEM_HH_ */
