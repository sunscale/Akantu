/**
 * @file   shape_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  shape functions for interface-enriched generalized FEM
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "shape_functions.hh"

#ifndef AKANTU_SHAPE_IGFEM_HH_
#define AKANTU_SHAPE_IGFEM_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <> class ShapeLagrange<_ek_igfem> : public ShapeFunctions {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ShapeLagrange(const Mesh & mesh, const ID & id = "shape_igfem",
                const MemoryID & memory_id = 0);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline void initShapeFunctions(const Array<Real> & nodes,
                                 const Matrix<Real> & integration_points,
                                 const Matrix<Real> & integration_points_1,
                                 const Matrix<Real> & integration_points_2,
                                 ElementType type,
                                 GhostType ghost_type);

  inline void
  interpolateEnrichmentsAllTypes(const Array<Real> & src, Array<Real> & dst,
                                 ElementType type,
                                 GhostType ghost_type) const;

  template <ElementType type>
  inline void precomputeShapesOnEnrichedNodes(const Array<Real> & nodes,
                                              GhostType ghost_type);

  template <ElementType type>
  void interpolateAtEnrichedNodes(const Array<Real> & src, Array<Real> & dst,
                                  GhostType ghost_type) const;

  /// pre compute all shapes on the element integration points from natural
  /// coordinates
  template <ElementType type>
  void precomputeShapesOnIntegrationPoints(const Array<Real> & nodes,
                                           GhostType ghost_type);

  /// pre compute all shape derivatives on the element integration points from
  /// natural coordinates
  template <ElementType type>
  void precomputeShapeDerivativesOnIntegrationPoints(const Array<Real> & nodes,
                                                     GhostType ghost_type);

  /// interpolate nodal values on the integration points
  template <ElementType type>
  void interpolateOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & uq, UInt nb_degree_of_freedom,
      GhostType ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const;

  /// interpolate on physical point
  template <ElementType type>
  void interpolate(const Vector<Real> & real_coords, UInt elem,
                   const Matrix<Real> & nodal_values,
                   Vector<Real> & interpolated,
                   GhostType ghost_type) const;

  /// compute the gradient of u on the integration points
  template <ElementType type>
  void gradientOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & nablauq, UInt nb_degree_of_freedom,
      GhostType ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const;

  /// multiply a field by shape functions  @f$ fts_{ij} = f_i * \varphi_j @f$
  template <ElementType type>
  void fieldTimesShapes(const Array<Real> & field,
                        Array<Real> & field_times_shapes,
                        GhostType ghost_type) const;

  /// find natural coords in parent element from real coords provided an element
  template <ElementType type>
  void inverseMap(const Vector<Real> & real_coords, UInt element,
                  Vector<Real> & natural_coords,
                  GhostType ghost_type = _not_ghost) const;

  /// find natural coords in sub-element from real coords provided an element
  template <ElementType type>
  void inverseMap(const Vector<Real> & real_coords, UInt element,
                  Vector<Real> & natural_coords, UInt sub_element,
                  GhostType ghost_type = _not_ghost) const;

  /// return true if the coordinates provided are inside the element, false
  /// otherwise
  template <ElementType type>
  bool contains(const Vector<Real> & real_coords, UInt elem,
                GhostType ghost_type) const;

  /// compute the shape on a provided point
  template <ElementType type>
  void computeShapes(const Vector<Real> & real_coords, UInt elem,
                     Vector<Real> & shapes, GhostType ghost_type) const;

  /// compute the shape derivatives on a provided point
  template <ElementType type>
  void computeShapeDerivatives(const Matrix<Real> & real_coords, UInt elem,
                               Tensor3<Real> & shapes,
                               GhostType ghost_type) const;

  /// interpolate a field on a given physical point
  template <ElementType type>
  void interpolateOnPhysicalPoint(const Vector<Real> & real_coords, UInt elem,
                                  const Array<Real> & field,
                                  Vector<Real> & interpolated,
                                  GhostType ghost_type) const;

  /// function to extract values at standard nodes and zero-out enriched values
  /// of a nodal field
  void extractValuesAtStandardNodes(const Array<Real> & nodal_values,
                                    Array<Real> & extracted_values,
                                    GhostType ghost_type) const;

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// compute the shape derivatives on integration points for a given element
  template <ElementType type>
  inline void
  computeShapeDerivativesOnCPointsByElement(const Matrix<Real> & node_coords,
                                            const Matrix<Real> & natural_coords,
                                            Tensor3<Real> & shapesd) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get a the shapes vector
  inline const Array<Real> &
  getShapes(ElementType el_type,
            GhostType ghost_type = _not_ghost) const;

  /// get a the shapes derivatives vector
  inline const Array<Real> &
  getShapesDerivatives(ElementType el_type,
                       GhostType ghost_type = _not_ghost) const;

  /// get a the shapes vector
  inline const Array<Real> &
  getShapesAtEnrichedNodes(ElementType el_type,
                           GhostType ghost_type = _not_ghost) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// shape functions for all elements
  ElementTypeMapArray<Real, InterpolationType> shapes;

  /// shape functions derivatives for all elements
  ElementTypeMapArray<Real, InterpolationType> shapes_derivatives;

  /// additional integration points for the IGFEM formulation
  ElementTypeMapArray<Real> igfem_integration_points;

  /// values of shape functions for all elements on the enriched nodes
  ElementTypeMapArray<Real, InterpolationType> shapes_at_enrichments;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "shape_igfem_inline_impl.hh"
/// standard output stream operator
// template <class ShapeFunction>
// inline std::ostream & operator <<(std::ostream & stream, const
// ShapeIGFEM<ShapeFunction> & _this)
// {
//   _this.printself(stream);
//   return stream;
// }

#endif /* AKANTU_SHAPE_IGFEM_HH_ */
