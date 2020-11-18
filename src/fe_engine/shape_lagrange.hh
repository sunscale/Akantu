/**
 * @file   shape_lagrange.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Feb 15 2011
 * @date last modification: Mon Jan 29 2018
 *
 * @brief  lagrangian shape functions class
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
#include "shape_lagrange_base.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SHAPE_LAGRANGE_HH_
#define AKANTU_SHAPE_LAGRANGE_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Shape> class ShapeCohesive;
class ShapeIGFEM;

template <ElementKind kind> class ShapeLagrange : public ShapeLagrangeBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ShapeLagrange(const Mesh & mesh, UInt spatial_dimension,
                const ID & id = "shape_lagrange",
                const MemoryID & memory_id = 0);
  ~ShapeLagrange() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialization function for structural elements not yet implemented
  inline void initShapeFunctions(const Array<Real> & nodes,
                                 const Matrix<Real> & integration_points,
                                 ElementType type,
                                 GhostType ghost_type);

  /// computes the shape functions derivatives for given interpolation points
  template <ElementType type>
  void computeShapeDerivativesOnIntegrationPoints(
      const Array<Real> & nodes, const Matrix<Real> & integration_points,
      Array<Real> & shape_derivatives, GhostType ghost_type,
      const Array<UInt> & filter_elements = empty_filter) const;

  void computeShapeDerivativesOnIntegrationPoints(
      const Array<Real> & nodes, const Matrix<Real> & integration_points,
      Array<Real> & shape_derivatives, ElementType type,
      GhostType ghost_type,
      const Array<UInt> & filter_elements) const override;

  /// pre compute all shapes on the element integration points from natural
  /// coordinates
  template <ElementType type>
  void precomputeShapesOnIntegrationPoints(const Array<Real> & nodes,
                                           GhostType ghost_type);

  /// pre compute all shape derivatives on the element integration points from
  /// natural coordinates
  template <ElementType type>
  void
  precomputeShapeDerivativesOnIntegrationPoints(const Array<Real> & nodes,
                                                GhostType ghost_type);

  /// interpolate nodal values on the integration points
  template <ElementType type>
  void interpolateOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & uq, UInt nb_degree_of_freedom,
      GhostType ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const;

  template <ElementType type>
  void interpolateOnIntegrationPoints(
      const Array<Real> & in_u, Array<Real> & out_uq, UInt nb_degree_of_freedom,
      const Array<Real> & shapes, GhostType ghost_type = _not_ghost,
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

  template <ElementType type>
  void computeBtD(const Array<Real> & Ds, Array<Real> & BtDs,
                  GhostType ghost_type,
                  const Array<UInt> & filter_elements) const;

  template <ElementType type>
  void computeBtDB(const Array<Real> & Ds, Array<Real> & BtDBs, UInt order_d,
                   GhostType ghost_type,
                   const Array<UInt> & filter_elements) const;

  /// multiply a field by shape functions  @f$ fts_{ij} = f_i * \varphi_j @f$
  template <ElementType type>
  void computeNtb(const Array<Real> & bs, Array<Real> & Ntbs,
                  GhostType ghost_type,
                  const Array<UInt> & filter_elements = empty_filter) const;

  /// find natural coords from real coords provided an element
  template <ElementType type>
  void inverseMap(const Vector<Real> & real_coords, UInt element,
                  Vector<Real> & natural_coords,
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

protected:
  /// compute the shape derivatives on integration points for a given element
  template <ElementType type>
  inline void
  computeShapeDerivativesOnCPointsByElement(const Matrix<Real> & node_coords,
                                            const Matrix<Real> & natural_coords,
                                            Tensor3<Real> & shapesd) const;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "shape_lagrange_inline_impl.hh"

#endif /* AKANTU_SHAPE_LAGRANGE_HH_ */
