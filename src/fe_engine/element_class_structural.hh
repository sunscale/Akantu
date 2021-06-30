/**
 * @file   element_class_structural.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Specialization of the element classes for structural elements
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_CLASS_STRUCTURAL_HH_
#define AKANTU_ELEMENT_CLASS_STRUCTURAL_HH_

namespace akantu {

/// Macro to generate the InterpolationProperty structures for different
/// interpolation types
#define AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(                  \
    itp_type, itp_geom_type, ndof, nb_stress, nb_dnds_cols)                    \
  template <> struct InterpolationProperty<itp_type> {                         \
    static const InterpolationKind kind{_itk_structural};                      \
    static const UInt nb_nodes_per_element{                                    \
        InterpolationProperty<itp_geom_type>::nb_nodes_per_element};           \
    static const InterpolationType itp_geometry_type{itp_geom_type};           \
    static const UInt natural_space_dimension{                                 \
        InterpolationProperty<itp_geom_type>::natural_space_dimension};        \
    static const UInt nb_degree_of_freedom{ndof};                              \
    static const UInt nb_stress_components{nb_stress};                         \
    static const UInt dnds_columns{nb_dnds_cols};                              \
  }

/* -------------------------------------------------------------------------- */
template <InterpolationType interpolation_type>
class InterpolationElement<interpolation_type, _itk_structural> {
public:
  using interpolation_property = InterpolationProperty<interpolation_type>;

  /// compute the shape values for a given set of points in natural coordinates
  static inline void computeShapes(const Matrix<Real> & natural_coord,
                                   const Matrix<Real> & real_coord,
                                   const Matrix<Real> & T, Tensor3<Real> & Ns) {
    for (UInt i = 0; i < natural_coord.cols(); ++i) {
      Matrix<Real> N_T = Ns(i);
      Matrix<Real> N(N_T.rows(), N_T.cols());
      computeShapes(natural_coord(i), real_coord, N);
      N_T.mul<false, false>(N, T);
    }
  }

  /// compute the shape values for a given point in natural coordinates
  static inline void computeShapes(const Vector<Real> & natural_coord,
                                   const Matrix<Real> & real_coord,
                                   Matrix<Real> & N);

  static inline void computeShapesMass(const Matrix<Real> & natural_coords,
                                       const Matrix<Real> & xs,
                                       const Matrix<Real> & T,
                                       Tensor3<Real> & Ns) {
    for (UInt i = 0; i < natural_coords.cols(); ++i) {
      Matrix<Real> N_T = Ns(i);
      Vector<Real> X = natural_coords(i);
      Matrix<Real> N(interpolation_property::nb_degree_of_freedom, N_T.cols());

      computeShapes(X, xs, N);
      N_T.mul<false, false>(N.block(0, 0, N_T.rows(), N_T.cols()), T);
    }
  }

  /// compute shape derivatives (input is dxds) for a set of points
  static inline void computeShapeDerivatives(const Tensor3<Real> & Js,
                                             const Tensor3<Real> & DNDSs,
                                             const Matrix<Real> & R,
                                             Tensor3<Real> & Bs) {
    for (UInt i = 0; i < Js.size(2); ++i) {
      Matrix<Real> J = Js(i);
      Matrix<Real> DNDS = DNDSs(i);
      Matrix<Real> DNDX(DNDS.rows(), DNDS.cols());
      auto inv_J = J.inverse();
      DNDX.mul<false, false>(inv_J, DNDS);
      Matrix<Real> B_R = Bs(i);
      Matrix<Real> B(B_R.rows(), B_R.cols());
      arrangeInVoigt(DNDX, B);
      B_R.mul<false, false>(B, R);
    }
  }

  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with variation of natural coordinates on a given set
   * of points in natural coordinates
   */
  static inline void computeDNDS(const Matrix<Real> & natural_coord,
                                 const Matrix<Real> & real_coord,
                                 Tensor3<Real> & dnds) {
    for (UInt i = 0; i < natural_coord.cols(); ++i) {
      Matrix<Real> dnds_t = dnds(i);
      computeDNDS(natural_coord(i), real_coord, dnds_t);
    }
  }

  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with
   * variation of natural coordinates on a given point in natural
   * coordinates
   */
  static inline void computeDNDS(const Vector<Real> & natural_coord,
                                 const Matrix<Real> & real_coord,
                                 Matrix<Real> & dnds);

  /**
   * arrange B in Voigt notation from DNDS
   */
  static inline void arrangeInVoigt(const Matrix<Real> & dnds,
                                    Matrix<Real> & B) {
    // Default implementation assumes dnds is already in Voigt notation
    B.deepCopy(dnds);
  }

public:
  static inline constexpr auto getNbNodesPerInterpolationElement() {
    return interpolation_property::nb_nodes_per_element;
  }

  static inline constexpr auto getShapeSize() {
    return interpolation_property::nb_nodes_per_element *
           interpolation_property::nb_degree_of_freedom *
           interpolation_property::nb_degree_of_freedom;
  }
  static inline constexpr auto getShapeIndependantSize() {
    return interpolation_property::nb_nodes_per_element *
           interpolation_property::nb_degree_of_freedom *
           interpolation_property::nb_stress_components;
  }
  static inline constexpr auto getShapeDerivativesSize() {
    return interpolation_property::nb_nodes_per_element *
           interpolation_property::nb_degree_of_freedom *
           interpolation_property::nb_stress_components;
  }
  static inline constexpr auto getNaturalSpaceDimension() {
    return interpolation_property::natural_space_dimension;
  }
  static inline constexpr auto getNbDegreeOfFreedom() {
    return interpolation_property::nb_degree_of_freedom;
  }
  static inline constexpr auto getNbStressComponents() {
    return interpolation_property::nb_stress_components;
  }
};

/// Macro to generate the element class structures for different structural
/// element types
/* -------------------------------------------------------------------------- */
#define AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(                       \
    elem_type, geom_type, interp_type, parent_el_type, elem_kind, sp,          \
    gauss_int_type, min_int_order)                                             \
  template <> struct ElementClassProperty<elem_type> {                         \
    static const GeometricalType geometrical_type{geom_type};                  \
    static const InterpolationType interpolation_type{interp_type};            \
    static const ElementType parent_element_type{parent_el_type};              \
    static const ElementKind element_kind{elem_kind};                          \
    static const UInt spatial_dimension{sp};                                   \
    static const GaussIntegrationType gauss_integration_type{gauss_int_type};  \
    static const UInt polynomial_degree{min_int_order};                        \
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
  using geometrical_element =
      GeometricalElement<ElementClassProperty<element_type>::geometrical_type>;
  using interpolation_element = InterpolationElement<
      ElementClassProperty<element_type>::interpolation_type>;
  using parent_element =
      ElementClass<ElementClassProperty<element_type>::parent_element_type>;

public:
  static inline void
  computeRotationMatrix(Matrix<Real> & /*R*/, const Matrix<Real> & /*X*/,
                        const Vector<Real> & /*extra_normal*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute jacobian (or integration variable change factor) for a given point
  static inline void computeJMat(const Vector<Real> & natural_coords,
                                 const Matrix<Real> & Xs, Matrix<Real> & J) {
    Matrix<Real> dnds(Xs.rows(), Xs.cols());
    parent_element::computeDNDS(natural_coords, dnds);
    J.mul<false, true>(dnds, Xs);
  }

  static inline void computeJMat(const Matrix<Real> & natural_coords,
                                 const Matrix<Real> & Xs, Tensor3<Real> & Js) {
    for (UInt i = 0; i < natural_coords.cols(); ++i) {
      // because non-const l-value reference does not bind to r-value
      Matrix<Real> J = Js(i);
      computeJMat(Vector<Real>(natural_coords(i)), Xs, J);
    }
  }

  static inline void computeJacobian(const Matrix<Real> & natural_coords,
                                     const Matrix<Real> & node_coords,
                                     Vector<Real> & jacobians) {
    using itp = typename interpolation_element::interpolation_property;
    Tensor3<Real> Js(itp::natural_space_dimension, itp::natural_space_dimension,
                     natural_coords.cols());
    computeJMat(natural_coords, node_coords, Js);
    for (UInt i = 0; i < natural_coords.cols(); ++i) {
      Matrix<Real> J = Js(i);
      jacobians(i) = J.det();
    }
  }

  static inline void computeRotation(const Matrix<Real> & node_coords,
                                     Matrix<Real> & rotation);

public:
  static AKANTU_GET_MACRO_NOT_CONST(Kind, _ek_structural, ElementKind);
  static AKANTU_GET_MACRO_NOT_CONST(P1ElementType, _not_defined, ElementType);
  static AKANTU_GET_MACRO_NOT_CONST(FacetType, _not_defined, ElementType);
  static constexpr auto getFacetType(__attribute__((unused)) UInt t = 0) {
    return _not_defined;
  }
  static constexpr AKANTU_GET_MACRO_NOT_CONST(
      SpatialDimension, ElementClassProperty<element_type>::spatial_dimension,
      UInt);
  static constexpr auto getFacetTypes() {
    return ElementClass<_not_defined>::getFacetTypes();
  }
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "element_class_hermite_inline_impl.hh"
/* keep order */
#include "element_class_bernoulli_beam_inline_impl.hh"
#include "element_class_kirchhoff_shell_inline_impl.hh"
/* -------------------------------------------------------------------------- */

#endif /* AKANTU_ELEMENT_CLASS_STRUCTURAL_HH_ */
