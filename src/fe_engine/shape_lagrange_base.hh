/**
 * @file   shape_lagrange_base.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 09 2017
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Base class for the shape lagrange
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "shape_functions.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SHAPE_LAGRANGE_BASE_HH_
#define AKANTU_SHAPE_LAGRANGE_BASE_HH_

namespace akantu {

class ShapeLagrangeBase : public ShapeFunctions {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ShapeLagrangeBase(const Mesh & mesh, UInt spatial_dimension,
                    ElementKind kind, const ID & id = "shape_lagrange",
                    const MemoryID & memory_id = 0);
  ~ShapeLagrangeBase() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// computes the shape functions for given interpolation points
  virtual void computeShapesOnIntegrationPoints(
      const Array<Real> & nodes, const Matrix<Real> & integration_points,
      Array<Real> & shapes, ElementType type,
      GhostType ghost_type,
      const Array<UInt> & filter_elements = empty_filter) const;

  /// computes the shape functions derivatives for given interpolation points
  virtual void computeShapeDerivativesOnIntegrationPoints(
      const Array<Real> & nodes, const Matrix<Real> & integration_points,
      Array<Real> & shape_derivatives, ElementType type,
      GhostType ghost_type,
      const Array<UInt> & filter_elements = empty_filter) const = 0;

  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  template <ElementType type>
  void computeShapesOnIntegrationPoints(
      const Array<Real> & nodes, const Matrix<Real> & integration_points,
      Array<Real> & shapes, GhostType ghost_type,
      const Array<UInt> & filter_elements = empty_filter) const;

public:
  void onElementsAdded(const Array<Element> & elements) override;
  void
  onElementsRemoved(const Array<Element> & elements,
                    const ElementTypeMapArray<UInt> & new_numbering) override;

protected:
  /// The kind to consider
  ElementKind _kind;
};

} // namespace akantu

#include "shape_lagrange_base_inline_impl.hh"

#endif /* AKANTU_SHAPE_LAGRANGE_BASE_HH_ */
