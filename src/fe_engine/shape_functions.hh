/**
 * @file   shape_functions.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  shape function class
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
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SHAPE_FUNCTIONS_HH_
#define AKANTU_SHAPE_FUNCTIONS_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
class ShapeFunctions  {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ShapeFunctions(const Mesh & mesh, UInt spatial_dimension,
                 const ID & id = "shape");
  virtual ~ShapeFunctions() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for (Int i = 0; i < indent; i++, space += AKANTU_INDENT) {
      ;
    }
    stream << space << "Shapes [" << std::endl;
    integration_points.printself(stream, indent + 1);
    // shapes.printself(stream, indent + 1);
    // shapes_derivatives.printself(stream, indent + 1);
    stream << space << "]" << std::endl;
  }

  /// set the integration points for a given element
  template <ElementType type>
  void setIntegrationPointsByType(const Matrix<Real> & integration_points,
                                  GhostType ghost_type);

  /// Build pre-computed matrices for interpolation of field form integration
  /// points at other given positions (interpolation_points)
  void initElementalFieldInterpolationFromIntegrationPoints(
      const ElementTypeMapArray<Real> & interpolation_points_coordinates,
      ElementTypeMapArray<Real> & interpolation_points_coordinates_matrices,
      ElementTypeMapArray<Real> & quad_points_coordinates_inv_matrices,
      const ElementTypeMapArray<Real> & quadrature_points_coordinates,
      const ElementTypeMapArray<UInt> * element_filter) const;

  /// Interpolate field at given position from given values of this field at
  /// integration points (field)
  /// using matrices precomputed with
  /// initElementalFieldInterplationFromIntegrationPoints
  void interpolateElementalFieldFromIntegrationPoints(
      const ElementTypeMapArray<Real> & field,
      const ElementTypeMapArray<Real> &
          interpolation_points_coordinates_matrices,
      const ElementTypeMapArray<Real> & quad_points_coordinates_inv_matrices,
      ElementTypeMapArray<Real> & result, GhostType ghost_type,
      const ElementTypeMapArray<UInt> * element_filter) const;

protected:
  /// interpolate nodal values stored by element on the integration points
  template <ElementType type>
  void interpolateElementalFieldOnIntegrationPoints(
      const Array<Real> & u_el, Array<Real> & uq, GhostType ghost_type,
      const Array<Real> & shapes,
      const Array<UInt> & filter_elements = empty_filter) const;

  /// gradient of nodal values stored by element on the control points
  template <ElementType type>
  void gradientElementalFieldOnIntegrationPoints(
      const Array<Real> & u_el, Array<Real> & out_nablauq,
      GhostType ghost_type, const Array<Real> & shapes_derivatives,
      const Array<UInt> & filter_elements) const;

protected:
  /// By element versions of non-templated eponym methods
  template <ElementType type>
  inline void interpolateElementalFieldFromIntegrationPoints(
      const Array<Real> & field,
      const Array<Real> & interpolation_points_coordinates_matrices,
      const Array<Real> & quad_points_coordinates_inv_matrices,
      ElementTypeMapArray<Real> & result, GhostType ghost_type,
      const Array<UInt> & element_filter) const;

  /// Interpolate field at given position from given values of this field at
  /// integration points (field)
  /// using matrices precomputed with
  /// initElementalFieldInterplationFromIntegrationPoints
  template <ElementType type>
  inline void initElementalFieldInterpolationFromIntegrationPoints(
      const Array<Real> & interpolation_points_coordinates,
      ElementTypeMapArray<Real> & interpolation_points_coordinates_matrices,
      ElementTypeMapArray<Real> & quad_points_coordinates_inv_matrices,
      const Array<Real> & quadrature_points_coordinates,
      GhostType ghost_type, const Array<UInt> & element_filter) const;

  /// build matrix for the interpolation of field form integration points
  template <ElementType type>
  inline void buildElementalFieldInterpolationMatrix(
      const Matrix<Real> & coordinates, Matrix<Real> & coordMatrix,
      UInt integration_order =
          ElementClassProperty<type>::polynomial_degree) const;

  /// build the so called interpolation matrix (first collumn is 1, then the
  /// other collumns are the traansposed coordinates)
  static inline void buildInterpolationMatrix(const Matrix<Real> & coordinates,
                                       Matrix<Real> & coordMatrix,
                                       UInt integration_order);

public:
  virtual void onElementsAdded(const Array<Element> & /*unused*/) {
    AKANTU_TO_IMPLEMENT();
  }
  virtual void onElementsRemoved(const Array<Element> & /*unused*/,
                                 const ElementTypeMapArray<UInt> & /*unused*/) {
    AKANTU_TO_IMPLEMENT();
  }
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the size of the shapes returned by the element class
  static inline UInt getShapeSize(ElementType type);

  /// get the size of the shapes derivatives returned by the element class
  static inline UInt getShapeDerivativesSize(ElementType type);

  inline const Matrix<Real> &
  getIntegrationPoints(ElementType type,
                       GhostType ghost_type) const {
    return integration_points(type, ghost_type);
  }

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

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// shape functions for all elements
  ElementTypeMapArray<Real, InterpolationType> shapes;

  /// shape functions derivatives for all elements
  ElementTypeMapArray<Real, InterpolationType> shapes_derivatives;

  /// associated mesh
  const Mesh & mesh;

  // spatial dimension of the elements to consider
  UInt _spatial_dimension;

  /// shape functions for all elements
  ElementTypeMap<Matrix<Real>> integration_points;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const ShapeFunctions & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu
#include "shape_functions_inline_impl.hh"

#endif /* AKANTU_SHAPE_FUNCTIONS_HH_ */
