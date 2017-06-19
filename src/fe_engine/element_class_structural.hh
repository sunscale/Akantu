/**
 * @file   element_class_structural.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Thu Oct 22 2015
 *
 * @brief  Specialization of the element classes for structural elements
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <InterpolationType interpolation_type>
class InterpolationElement<interpolation_type, _itk_structural> {
public:
  typedef InterpolationPorperty<interpolation_type> interpolation_property;

  /// compute the shape values for a given set of points in natural coordinates
  static inline void computeShapes(const Matrix<Real> & natural_coord,
                                   Matrix<Real> & N,
                                   const Matrix<Real> & real_nodal_coord,
                                   UInt n = 0) {
    UInt nb_points = natural_coord.cols();
    for (UInt p = 0; p < nb_points; ++p) {
      Vector<Real> Np = N(p);
      computeShapes(natural_coord(p), Np, real_nodal_coord, n);
    }
  }

  /// compute the shape values for a given point in natural coordinates
  static inline void computeShapes(const Vector<Real> & natural_coord,
                                   Vector<Real> & N,
                                   const Matrix<Real> & real_nodal_coord,
                                   UInt n = 0);
  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with variation of natural coordinates on a given set
   * of points in natural coordinates
   */
  static inline void computeDNDS(const Matrix<Real> & natural_coord,
                                 Tensor3<Real> & dnds,
                                 const Matrix<Real> & real_nodal_coord,
                                 UInt n = 0) {
    for (UInt i = 0; i < natural_coord.cols(); ++i) {
      Matrix<Real> dnds_t = dnds(i);
      computeDNDS(natural_coord(i), dnds_t, real_nodal_coord, n);
    }
  }

  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with
   * variation of natural coordinates on a given point in natural
   * coordinates
   */
  static inline void computeDNDS(const Vector<Real> & natural_coord,
                                 Matrix<Real> & dnds,
                                 const Matrix<Real> & real_nodal_coord,
                                 UInt n = 0);

public:
  static AKANTU_GET_MACRO_NOT_CONST(NbShapeFunctions, nb_shape_functions, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NbShapeDerivatives, nb_shape_derivatives,
                                    UInt);
  static AKANTU_GET_MACRO_NOT_CONST(
      ShapeSize, interpolation_property::nb_nodes_per_element, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(
      ShapeDerivativesSize, (interpolation_property::nb_nodes_per_element *
                             interpolation_property::natural_space_dimension),
      UInt);
  static AKANTU_GET_MACRO_NOT_CONST(
      NaturalSpaceDimension, interpolation_property::natural_space_dimension,
      UInt);

protected:
  /// nb shape functions
  static const UInt nb_shape_functions;
  static const UInt nb_shape_derivatives;
};

/// Macro to generate the element class structures for different structural
/// element types
/* -------------------------------------------------------------------------- */
#define AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(                       \
    elem_type, geom_type, interp_type, parent_el_type, elem_kind, sp,          \
    gauss_int_type, min_int_order)                                             \
  template <> struct ElementClassProperty<elem_type> {                         \
    static const GeometricalType geometrical_type = geom_type;                 \
    static const InterpolationType interpolation_type = interp_type;           \
    static const ElementType parent_element_type = parent_el_type;             \
    static const ElementKind element_kind = elem_kind;                         \
    static const UInt spatial_dimension = sp;                                  \
    static const GaussIntegrationType gauss_integration_type = gauss_int_type; \
    static const UInt polynomial_degree = min_int_order;                       \
  }

/* -------------------------------------------------------------------------- */
/* ElementClass for structural elements                                       */
/* -------------------------------------------------------------------------- */
template <ElementType element_type>
class ElementClass<element_type, _ek_structural>
    : public GeometricalElement<
          ElementClassProperty<element_type>::geometrical_type>,
      public InterpolationElement<
          ElementClassProperty<element_type>::interpolation_type> {
protected:
  typedef GeometricalElement<
      ElementClassProperty<element_type>::geometrical_type>
      geometrical_element;
  typedef InterpolationElement<
      ElementClassProperty<element_type>::interpolation_type>
      interpolation_element;
  typedef ElementClass<ElementClassProperty<element_type>::parent_element_type>
      parent_element;

public:
  /// compute shape derivatives (input is dxds) for a set of points
  static inline void
  computeShapeDerivatives(const Matrix<Real> & natural_coord,
                          Tensor3<Real> & shape_deriv,
                          const Matrix<Real> & real_nodal_coord, UInt n = 0) {
    UInt nb_points = natural_coord.cols();

    if (element_type == _kirchhoff_shell) {
      /// TO BE CONTINUED and moved in a _tmpl.hh

      UInt spatial_dimension = real_nodal_coord.cols();
      UInt nb_nodes = real_nodal_coord.rows();

      const UInt projected_dim = natural_coord.rows();
      Matrix<Real> rotation_matrix(real_nodal_coord);
      Matrix<Real> rotated_nodal_coord(real_nodal_coord);
      Matrix<Real> projected_nodal_coord(natural_coord);
      /* --------------------------------------------------------------------------
       */
      /* --------------------------------------------------------------------------
       */

      Matrix<Real> Pe(real_nodal_coord);
      Matrix<Real> Pg(real_nodal_coord);
      Matrix<Real> inv_Pg(real_nodal_coord);

      /// compute matrix Pe
      Pe.eye();

      /// compute matrix Pg
      Vector<Real> Pg_col_1(spatial_dimension);
      Pg_col_1(0) = real_nodal_coord(0, 1) - real_nodal_coord(0, 0);
      Pg_col_1(1) = real_nodal_coord(1, 1) - real_nodal_coord(1, 0);
      Pg_col_1(2) = real_nodal_coord(2, 1) - real_nodal_coord(2, 0);

      Vector<Real> Pg_col_2(spatial_dimension);
      Pg_col_2(0) = real_nodal_coord(0, 2) - real_nodal_coord(0, 0);
      Pg_col_2(1) = real_nodal_coord(1, 2) - real_nodal_coord(1, 0);
      Pg_col_2(2) = real_nodal_coord(2, 2) - real_nodal_coord(2, 0);

      Vector<Real> Pg_col_3(spatial_dimension);
      Pg_col_3.crossProduct(Pg_col_1, Pg_col_2);

      for (UInt i = 0; i < nb_points; ++i) {
        Pg(i, 0) = Pg_col_1(i);
        Pg(i, 1) = Pg_col_2(i);
        Pg(i, 2) = Pg_col_3(i);
      }

      /// compute inverse of Pg
      inv_Pg.inverse(Pg);

      /// compute rotation matrix
      // rotation_matrix=Pe*inv_Pg;
      rotation_matrix.eye();
      /* --------------------------------------------------------------------------
       */
      /* --------------------------------------------------------------------------
       */

      rotated_nodal_coord.mul<false, false>(rotation_matrix, real_nodal_coord);

      for (UInt i = 0; i < projected_dim; ++i) {
        for (UInt j = 0; j < nb_points; ++j) {
          projected_nodal_coord(i, j) = rotated_nodal_coord(i, j);
        }
      }

      Tensor3<Real> dnds(projected_dim, nb_nodes, natural_coord.cols());
      Tensor3<Real> J(projected_dim, projected_dim, natural_coord.cols());

      parent_element::computeDNDS(natural_coord, dnds);

      parent_element::computeJMat(dnds, projected_nodal_coord, J);

      for (UInt p = 0; p < nb_points; ++p) {

        Matrix<Real> shape_deriv_p = shape_deriv(p);

        interpolation_element::computeDNDS(natural_coord(p), shape_deriv_p,
                                           projected_nodal_coord, n);

        Matrix<Real> dNdS = shape_deriv_p;
        Matrix<Real> inv_J(projected_dim, projected_dim);
        inv_J.inverse(J(p));
        shape_deriv_p.mul<false, false>(inv_J, dNdS);
      }
    } else {

      for (UInt p = 0; p < nb_points; ++p) {
        Matrix<Real> shape_deriv_p = shape_deriv(p);
        interpolation_element::computeDNDS(natural_coord(p), shape_deriv_p,
                                           real_nodal_coord, n);
      }
    }
  }

  /// compute jacobian (or integration variable change factor) for a given point
  static inline void computeJacobian(const Matrix<Real> & natural_coords,
                                     const Matrix<Real> & nodal_coords,
                                     Vector<Real> & jacobians) {
    parent_element::computeJacobian(natural_coords, nodal_coords, jacobians);
  }

public:
  static AKANTU_GET_MACRO_NOT_CONST(Kind, _ek_structural, ElementKind);
  static AKANTU_GET_MACRO_NOT_CONST(P1ElementType, _not_defined, ElementType);
  static AKANTU_GET_MACRO_NOT_CONST(FacetType, _not_defined, ElementType);
  static ElementType getFacetType(__attribute__((unused)) UInt t = 0) {
    return _not_defined;
  }
  static AKANTU_GET_MACRO_NOT_CONST(
      SpatialDimension, ElementClassProperty<element_type>::spatial_dimension,
      UInt);
  static ElementType * getFacetTypeInternal() { return NULL; }
};

#include "element_classes/element_class_bernoulli_beam_inline_impl.cc"
#include "element_classes/element_class_kirchhoff_shell_inline_impl.cc"

} // akantu
