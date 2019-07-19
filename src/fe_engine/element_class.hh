/**
 * @file   element_class.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Declaration of the ElementClass main class and the
 * Integration and Interpolation elements
 *
 * @section LICENSE
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
#include "aka_common.hh"
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_CLASS_HH__
#define __AKANTU_ELEMENT_CLASS_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
/// default element class structure
template <ElementType element_type> struct ElementClassProperty {
  static const GeometricalType geometrical_type{_gt_not_defined};
  static const InterpolationType interpolation_type{_itp_not_defined};
  static const ElementKind element_kind{_ek_regular};
  static const UInt spatial_dimension{0};
  static const GaussIntegrationType gauss_integration_type{_git_not_defined};
  static const UInt polynomial_degree{0};
};

/// Macro to generate the element class structures for different element types
#define AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(elem_type, geom_type,             \
                                             interp_type, elem_kind, sp,       \
                                             gauss_int_type, min_int_order)    \
  template <> struct ElementClassProperty<elem_type> {                         \
    static const GeometricalType geometrical_type{geom_type};                  \
    static const InterpolationType interpolation_type{interp_type};            \
    static const ElementKind element_kind{elem_kind};                          \
    static const UInt spatial_dimension{sp};                                   \
    static const GaussIntegrationType gauss_integration_type{gauss_int_type};  \
    static const UInt polynomial_degree{min_int_order};                        \
  }

/* -------------------------------------------------------------------------- */
/* Geometry                                                                   */
/* -------------------------------------------------------------------------- */
/// Default GeometricalShape structure
template <GeometricalType geometrical_type> struct GeometricalShape {
  static const GeometricalShapeType shape{_gst_point};
};

/// Templated GeometricalShape with function contains
template <GeometricalShapeType shape> struct GeometricalShapeContains {
  /// Check if the point (vector in 2 and 3D) at natural coordinate coor
  template <class vector_type>
  static inline bool contains(const vector_type & coord);
};

/// Macro to generate the GeometricalShape structures for different geometrical
/// types
#define AKANTU_DEFINE_SHAPE(geom_type, geom_shape)                             \
  template <> struct GeometricalShape<geom_type> {                             \
    static const GeometricalShapeType shape{geom_shape};                       \
  }

AKANTU_DEFINE_SHAPE(_gt_hexahedron_20, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_hexahedron_8, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_pentahedron_15, _gst_prism);
AKANTU_DEFINE_SHAPE(_gt_pentahedron_6, _gst_prism);
AKANTU_DEFINE_SHAPE(_gt_point, _gst_point);
AKANTU_DEFINE_SHAPE(_gt_quadrangle_4, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_quadrangle_8, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_segment_2, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_segment_3, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_tetrahedron_10, _gst_triangle);
AKANTU_DEFINE_SHAPE(_gt_tetrahedron_4, _gst_triangle);
AKANTU_DEFINE_SHAPE(_gt_triangle_3, _gst_triangle);
AKANTU_DEFINE_SHAPE(_gt_triangle_6, _gst_triangle);

/* -------------------------------------------------------------------------- */
template <GeometricalType geometrical_type>
struct GeometricalElementProperty {};

template <ElementType element_type>
struct ElementClassExtraGeometryProperties {};

/* -------------------------------------------------------------------------- */
/// Templated GeometricalElement with function getInradius
template <GeometricalType geometrical_type,
          GeometricalShapeType shape =
              GeometricalShape<geometrical_type>::shape>
class GeometricalElement {
  using geometrical_property = GeometricalElementProperty<geometrical_type>;

public:
  /// compute the in-radius: \todo should be renamed for characteristic length
  static inline Real getInradius(__attribute__((unused))
                                 const Matrix<Real> & coord) {
    AKANTU_TO_IMPLEMENT();
  }

  /// true if the natural coordinates are in the element
  template <class vector_type>
  static inline bool contains(const vector_type & coord);

public:
  static AKANTU_GET_MACRO_NOT_CONST(SpatialDimension,
                                    geometrical_property::spatial_dimension,
                                    UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NbNodesPerElement,
                                    geometrical_property::nb_nodes_per_element,
                                    UInt);
  static inline constexpr auto getNbFacetTypes() {
    return geometrical_property::nb_facet_types;
  };
  static inline UInt getNbFacetsPerElement(UInt t);
  static inline UInt getNbFacetsPerElement();
  static inline constexpr auto getFacetLocalConnectivityPerElement(UInt t = 0);
};

/* -------------------------------------------------------------------------- */
/* Interpolation                                                              */
/* -------------------------------------------------------------------------- */
/// default InterpolationProperty structure
template <InterpolationType interpolation_type> struct InterpolationProperty {};

/// Macro to generate the InterpolationProperty structures for different
/// interpolation types
#define AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(itp_type, itp_kind,          \
                                                  nb_nodes, ndim)              \
  template <> struct InterpolationProperty<itp_type> {                         \
    static constexpr InterpolationKind kind{itp_kind};                         \
    static constexpr UInt nb_nodes_per_element{nb_nodes};                      \
    static constexpr UInt natural_space_dimension{ndim};                       \
  }

/* -------------------------------------------------------------------------- */
/// Generic (templated by the enum InterpolationType which specifies the order
/// and the dimension of the interpolation) class handling the elemental
/// interpolation
template <InterpolationType interpolation_type,
          InterpolationKind kind =
              InterpolationProperty<interpolation_type>::kind>
class InterpolationElement {
public:
  using interpolation_property = InterpolationProperty<interpolation_type>;

  /// compute the shape values for a given set of points in natural coordinates
  static inline void computeShapes(const Matrix<Real> & natural_coord,
                                   Matrix<Real> & N);

  /// compute the shape values for a given point in natural coordinates
  template <class vector_type>
  static inline void computeShapes(const vector_type &, vector_type &) {
    AKANTU_TO_IMPLEMENT();
  }

  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with variation of natural coordinates on a given set
   * of points in natural coordinates
   */
  static inline void computeDNDS(const Matrix<Real> & natural_coord,
                                 Tensor3<Real> & dnds);
  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with
   * variation of natural coordinates on a given point in natural
   * coordinates
   */
  template <class vector_type, class matrix_type>
  static inline void computeDNDS(const vector_type &, matrix_type &) {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute jacobian (or integration variable change factor) for a given point
  /// in the case of spatial_dimension != natural_space_dimension
  static inline void computeSpecialJacobian(const Matrix<Real> &, Real &) {
    AKANTU_TO_IMPLEMENT();
  }

  /// interpolate a field given (arbitrary) natural coordinates
  static inline void
  interpolateOnNaturalCoordinates(const Vector<Real> & natural_coords,
                                  const Matrix<Real> & nodal_values,
                                  Vector<Real> & interpolated);

  /// interpolate a field given the shape functions on the interpolation point
  static inline void interpolate(const Matrix<Real> & nodal_values,
                                 const Vector<Real> & shapes,
                                 Vector<Real> & interpolated);

  /// interpolate a field given the shape functions on the interpolations points
  static inline void interpolate(const Matrix<Real> & nodal_values,
                                 const Matrix<Real> & shapes,
                                 Matrix<Real> & interpolated);

  /// compute the gradient of a given field on the given natural coordinates
  static inline void
  gradientOnNaturalCoordinates(const Vector<Real> & natural_coords,
                               const Matrix<Real> & f, Matrix<Real> & gradient);

public:
  static AKANTU_GET_MACRO_NOT_CONST(
      ShapeSize,
      InterpolationProperty<interpolation_type>::nb_nodes_per_element, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(
      ShapeDerivativesSize,
      (InterpolationProperty<interpolation_type>::nb_nodes_per_element *
       InterpolationProperty<interpolation_type>::natural_space_dimension),
      UInt);
  static AKANTU_GET_MACRO_NOT_CONST(
      NaturalSpaceDimension,
      InterpolationProperty<interpolation_type>::natural_space_dimension, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(
      NbNodesPerInterpolationElement,
      InterpolationProperty<interpolation_type>::nb_nodes_per_element, UInt);
};

/* -------------------------------------------------------------------------- */
/* Integration                                                                */
/* -------------------------------------------------------------------------- */
template <GaussIntegrationType git_class, UInt nb_points>
struct GaussIntegrationTypeData {
  /// quadrature points in natural coordinates
  static Real quad_positions[];
  /// weights for the Gauss integration
  static Real quad_weights[];
};

template <ElementType type,
          UInt n = ElementClassProperty<type>::polynomial_degree>
class GaussIntegrationElement {
public:
  static UInt getNbQuadraturePoints();
  static const Matrix<Real> getQuadraturePoints();
  static const Vector<Real> getWeights();
};

/* -------------------------------------------------------------------------- */
/* ElementClass                                                               */
/* -------------------------------------------------------------------------- */
template <ElementType element_type,
          ElementKind element_kind =
              ElementClassProperty<element_type>::element_kind>
class ElementClass
    : public GeometricalElement<
          ElementClassProperty<element_type>::geometrical_type>,
      public InterpolationElement<
          ElementClassProperty<element_type>::interpolation_type> {
protected:
  using geometrical_element =
      GeometricalElement<ElementClassProperty<element_type>::geometrical_type>;
  using interpolation_element = InterpolationElement<
      ElementClassProperty<element_type>::interpolation_type>;

  using element_property = ElementClassProperty<element_type>;
  using interpolation_property =
      typename interpolation_element::interpolation_property;

public:
  /**
   * compute @f$ J = \frac{\partial x_j}{\partial s_i} @f$ the variation of real
   * coordinates along with variation of natural coordinates on a given point in
   * natural coordinates
   */
  static inline void computeJMat(const Matrix<Real> & dnds,
                                 const Matrix<Real> & node_coords,
                                 Matrix<Real> & J);

  /**
   * compute the Jacobian matrix by computing the variation of real coordinates
   * along with variation of natural coordinates on a given set of points in
   * natural coordinates
   */
  static inline void computeJMat(const Tensor3<Real> & dnds,
                                 const Matrix<Real> & node_coords,
                                 Tensor3<Real> & J);

  /// compute the jacobians of a serie of natural coordinates
  static inline void computeJacobian(const Matrix<Real> & natural_coords,
                                     const Matrix<Real> & node_coords,
                                     Vector<Real> & jacobians);

  /// compute jacobian (or integration variable change factor) for a set of
  /// points
  static inline void computeJacobian(const Tensor3<Real> & J,
                                     Vector<Real> & jacobians);

  /// compute jacobian (or integration variable change factor) for a given point
  static inline void computeJacobian(const Matrix<Real> & J, Real & jacobians);

  /// compute shape derivatives (input is dxds) for a set of points
  static inline void computeShapeDerivatives(const Tensor3<Real> & J,
                                             const Tensor3<Real> & dnds,
                                             Tensor3<Real> & shape_deriv);

  /// compute shape derivatives (input is dxds) for a given point
  static inline void computeShapeDerivatives(const Matrix<Real> & J,
                                             const Matrix<Real> & dnds,
                                             Matrix<Real> & shape_deriv);

  /// compute the normal of a surface defined by the function f
  static inline void
  computeNormalsOnNaturalCoordinates(const Matrix<Real> & coord,
                                     Matrix<Real> & f, Matrix<Real> & normals);

  /// get natural coordinates from real coordinates
  static inline void inverseMap(const Vector<Real> & real_coords,
                                const Matrix<Real> & node_coords,
                                Vector<Real> & natural_coords,
                                Real tolerance = 1e-10);

  /// get natural coordinates from real coordinates
  static inline void inverseMap(const Matrix<Real> & real_coords,
                                const Matrix<Real> & node_coords,
                                Matrix<Real> & natural_coords,
                                Real tolerance = 1e-10);

public:
  static AKANTU_GET_MACRO_NOT_CONST(Kind, element_kind, ElementKind);
  static constexpr AKANTU_GET_MACRO_NOT_CONST(
      SpatialDimension, ElementClassProperty<element_type>::spatial_dimension,
      UInt);

  using element_class_extra_geom_property =
      ElementClassExtraGeometryProperties<element_type>;

  static constexpr auto getP1ElementType() {
    return element_class_extra_geom_property::p1_type;
  }
  static constexpr auto getFacetType(UInt t = 0) {
    return element_class_extra_geom_property::facet_type[t];
  }
  static constexpr auto getFacetTypes();
};

/* -------------------------------------------------------------------------- */
} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "geometrical_element_property.hh"
#include "interpolation_element_tmpl.hh"
/* -------------------------------------------------------------------------- */
#include "element_class_tmpl.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
#include "element_class_hexahedron_8_inline_impl.cc"
#include "element_class_pentahedron_6_inline_impl.cc"
/* keep order */
#include "element_class_hexahedron_20_inline_impl.cc"
#include "element_class_pentahedron_15_inline_impl.cc"
#include "element_class_point_1_inline_impl.cc"
#include "element_class_quadrangle_4_inline_impl.cc"
#include "element_class_quadrangle_8_inline_impl.cc"
#include "element_class_segment_2_inline_impl.cc"
#include "element_class_segment_3_inline_impl.cc"
#include "element_class_tetrahedron_10_inline_impl.cc"
#include "element_class_tetrahedron_4_inline_impl.cc"
#include "element_class_triangle_3_inline_impl.cc"
#include "element_class_triangle_6_inline_impl.cc"
} // namespace akantu

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_STRUCTURAL_MECHANICS)
#include "element_class_structural.hh"
#endif

#if defined(AKANTU_COHESIVE_ELEMENT)
#include "cohesive_element.hh"
#endif

#if defined(AKANTU_IGFEM)
#include "element_class_igfem.hh"
#endif

#endif /* __AKANTU_ELEMENT_CLASS_HH__ */
