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
/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_CLASS_STRUCTURAL_HH__
#define __AKANTU_ELEMENT_CLASS_STRUCTURAL_HH__

namespace akantu {

/// Macro to generate the InterpolationProperty structures for different
/// interpolation types
#define AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(                  \
    itp_type, itp_geom_type, ndof, nb_stress)                                  \
  template <> struct InterpolationProperty<itp_type> {                         \
    static const InterpolationKind kind{_itk_structural};                      \
    static const UInt nb_nodes_per_element{                                    \
        InterpolationProperty<itp_geom_type>::nb_nodes_per_element};           \
    static const InterpolationType itp_geometry_type{itp_geom_type};           \
    static const UInt natural_space_dimension{                                 \
        InterpolationProperty<itp_geom_type>::natural_space_dimension};        \
    static const UInt nb_degree_of_freedom{ndof};                              \
    static const UInt nb_stress_components{nb_stress};                         \
  }

/* -------------------------------------------------------------------------- */
template <InterpolationType interpolation_type>
class InterpolationElement<interpolation_type, _itk_structural> {
public:
  using interpolation_property = InterpolationProperty<interpolation_type>;

  /// compute the shape values for a given set of points in natural coordinates
  static inline void computeShapes(const Matrix<Real> & natural_coord,
                                   Tensor3<Real> & N) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// compute the shape values for a given point in natural coordinates
  static inline void computeShapes(const Vector<Real> & natural_coord,
                                   Matrix<Real> & N);
  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with variation of natural coordinates on a given set
   * of points in natural coordinates
   */
  static inline void computeDNDS(const Matrix<Real> & natural_coord,
                                 Tensor3<Real> & dnds) {
    for (UInt i = 0; i < natural_coord.cols(); ++i) {
      Matrix<Real> dnds_t = dnds(i);
      computeDNDS(natural_coord(i), dnds_t);
    }
  }

  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with
   * variation of natural coordinates on a given point in natural
   * coordinates
   */
  static inline void computeDNDS(const Vector<Real> & natural_coord,
                                 Matrix<Real> & dnds);

public:
  static AKANTU_GET_MACRO_NOT_CONST(
      NbNodesPerInterpolationElement,
      interpolation_property::nb_nodes_per_element, UInt);

  static AKANTU_GET_MACRO_NOT_CONST(
      ShapeSize,
      (interpolation_property::nb_nodes_per_element *
       interpolation_property::nb_degree_of_freedom *
       interpolation_property::nb_degree_of_freedom),
      UInt);
  static AKANTU_GET_MACRO_NOT_CONST(
      ShapeDerivativesSize,
      (interpolation_property::nb_nodes_per_element *
       interpolation_property::nb_degree_of_freedom *
       interpolation_property::nb_stress_components),
      UInt);
  static AKANTU_GET_MACRO_NOT_CONST(
      NaturalSpaceDimension, interpolation_property::natural_space_dimension,
      UInt);
  static AKANTU_GET_MACRO_NOT_CONST(
      NbDegreeOfFreedom, interpolation_property::nb_degree_of_freedom, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(
      VoigtSize, interpolation_property::voigt_size, UInt);

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
  /// compute shape derivatives (input is dxds) for a set of points
  static inline void
  computeShapeDerivatives(const Matrix<Real> & natural_coord,
                          Tensor3<Real> & shape_deriv,
                          const Matrix<Real> & real_nodal_coord) {
    UInt nb_points = natural_coord.cols();

    // compute dnds
    Tensor3<Real> dnds(real_nodal_coord.rows(), real_nodal_coord.cols(),
                       natural_coord.cols());
    ElementClass<element_type>::computeDNDS(natural_coord, dnds);
    // compute jacobian
    Tensor3<Real> J(real_nodal_coord.rows(), natural_coord.rows(),
                    natural_coord.cols());
    ElementClass<element_type>::computeJMat(dnds, real_nodal_coord, J);

    // compute dndx
    ElementClass<element_type>::computeShapeDerivatives(J, dnds, shape_deriv);

  }

  /// compute jacobian (or integration variable change factor) for a given point
  static inline void computeJMat(const Tensor3<Real> & natural_coords,
                                 const Matrix<Real> & nodal_coords,
                                 Tensor3<Real> & jacobian) {
    AKANTU_DEBUG_TO_IMPLEMENT();
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

} // namespace akantu

#include "element_classes/element_class_hermite_inline_impl.cc"

#include "element_classes/element_class_bernoulli_beam_inline_impl.cc"
//#include "element_classes/element_class_kirchhoff_shell_inline_impl.cc"

#endif /* __AKANTU_ELEMENT_CLASS_STRUCTURAL_HH__ */
