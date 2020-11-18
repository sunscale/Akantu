/**
 * @file   element_class_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Specialization for interface-enriched finite elements
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_CLASS_IGFEM_HH_
#define AKANTU_ELEMENT_CLASS_IGFEM_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <InterpolationType interpolation_type>
class InterpolationElement<interpolation_type, _itk_igfem> {
public:
  using interpolation_property = InterpolationProperty<interpolation_type>;
  /* ------------------------------------------------------------------------ */
  /* Member functions                                                         */
  /* ------------------------------------------------------------------------ */
public:
  static void assembleShapes(const Vector<Real> & parent_interpolation,
                             const Vector<Real> & sub_interpolation,
                             Vector<Real> & interpolation,
                             UInt sub_element = 0) {
    /// N1, N2, N3 of parent triangle
    UInt nb_nodes_parent = InterpolationElement<
        interpolation_property::parent_interpolation_type>::getShapeSize();

    for (UInt i = 0; i < nb_nodes_parent; ++i) {
      interpolation(i) = parent_interpolation(i);
    }
    /// add the enrichment
    UInt * enriched_node = enrichments[sub_element];
    for (UInt e = 0; e < nb_enrichments; ++e) {
      interpolation(nb_nodes_parent + e) = sub_interpolation(enriched_node[e]);
    }
  }
  static void
  assembleShapeDerivatives(const Matrix<Real> & parent_interpolation,
                           const Matrix<Real> & sub_interpolation,
                           Matrix<Real> & interpolation, UInt sub_element = 0) {

    /// N1, N2, N3 of parent triangle
    UInt nb_nodes_parent = InterpolationElement<
        interpolation_property::parent_interpolation_type>::getShapeSize();

    for (UInt i = 0; i < nb_nodes_parent; ++i) {
      Vector<Real> ip(interpolation(i));
      ip = parent_interpolation(i);
    }
    /// add the enrichment
    UInt * enriched_node = enrichments[sub_element];
    for (UInt e = 0; e < nb_enrichments; ++e) {
      Vector<Real> ip(interpolation(nb_nodes_parent + e));
      ip = sub_interpolation(enriched_node[e]);
    }
  }

  static void interpolate(const Matrix<Real> & nodal_values,
                          const Vector<Real> & shapes,
                          Vector<Real> & interpolated) {
    Matrix<Real> interpm(interpolated.storage(), nodal_values.rows(), 1);
    Matrix<Real> shapesm(shapes.storage(),
                         interpolation_property::nb_nodes_per_element, 1);
    interpm.mul<false, false>(nodal_values, shapesm);
  }

public:
  static AKANTU_GET_MACRO_NOT_CONST(
      ShapeSize, interpolation_property::nb_nodes_per_element, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(
      ShapeDerivativesSize,
      (interpolation_property::nb_nodes_per_element *
       interpolation_property::natural_space_dimension),
      UInt);
  static AKANTU_GET_MACRO_NOT_CONST(
      NaturalSpaceDimension, interpolation_property::natural_space_dimension,
      UInt);
  static AKANTU_GET_MACRO_NOT_CONST(
      NbNodesPerInterpolationElement,
      interpolation_property::nb_nodes_per_element, UInt);
  static AKANTU_GET_MACRO_NOT_CONST(NbSubElements, interpolation_property::nb_sub_elements, UInt);
  static UInt * getSubElementConnectivity(UInt t = 0) {
    return &(interpolation_property::sub_element_connectivity[t]);
  };
  static UInt getNbEnrichments() { return interpolation_property::nb_enrichments; };
  static UInt * getSubElementEnrichments(UInt t = 0) { return &(interpolation_property::enrichments[t]); };

protected:
  /// storage of the subelement local connectivity

  static UInt sub_element_connectivity_vect[];
  /// local connectivity of subelements
  static UInt * sub_element_connectivity[];
  /// nb of subelements
  static UInt nb_sub_elements;
  /// storage of enrichments
  static UInt enrichment_vect[];
  static UInt * enrichments[];
  static UInt nb_enrichments;
};

} // namespace akantu

#include "interpolation_element_igfem_tmpl.hh"

namespace akantu {
/* -------------------------------------------------------------------------- */
#define AKANTU_DEFINE_IGFEM_ELEMENT_CLASS_PROPERTY(                            \
    elem_type, geom_type, interp_type, parent_el_type, sub_el_type_1,          \
    sub_el_type_2, elem_kind, sp, min_int_order)                               \
  template <> struct ElementClassProperty<elem_type> {                         \
    static const GeometricalType geometrical_type{geom_type};                  \
    static const InterpolationType interpolation_type{interp_type};            \
    static const ElementType parent_element_type{parent_el_type};              \
    static const ElementType sub_element_type_1{sub_el_type_1};                \
    static const ElementType sub_element_type_2{sub_el_type_2};                \
    static const ElementKind element_kind{elem_kind};                          \
    static const UInt spatial_dimension{sp};                                   \
    static const UInt minimal_integration_order{min_int_order};                \
  }

/* -------------------------------------------------------------------------- */
template <ElementType element_type>
class ElementClass<element_type, _ek_igfem>
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
  using element_property = ElementClassProperty<element_type>;
  using interpolation_property =
      typename interpolation_element::interpolation_property;

  /* ------------------------------------------------------------------------ */
  /* Member functions                                                         */
  /* ------------------------------------------------------------------------ */
public:
  static void getSubElementCoords(const Matrix<Real> & element_coords,
                                  Matrix<Real> & sub_coords,
                                  const UInt sub_element) {
    /// get the sub_element_type
    ///   constexrp ElementType sub_el_type = getSubElementType(sub_element);
    UInt nb_nodes_sub_el = 0;
    switch (sub_element) {
    case 0:
      nb_nodes_sub_el =
          ElementClass<ElementClassProperty<element_type>::sub_element_type_1>::
              getNbNodesPerInterpolationElement();
      break;
    case 1:
      nb_nodes_sub_el =
          ElementClass<ElementClassProperty<element_type>::sub_element_type_2>::
              getNbNodesPerInterpolationElement();
      break;
    }

    for (UInt i = 0; i < nb_nodes_sub_el; ++i) {
      UInt lc = InterpolationElement<
          ElementClassProperty<element_type>::interpolation_type>::
          sub_element_connectivity[sub_element][i];
      Vector<Real> sub_c(sub_coords(i));
      sub_c = element_coords(lc);
    }
  }

  static void getParentCoords(const Matrix<Real> & element_coords,
                              Matrix<Real> & parent_coords) {
    const ElementType parent_type =
        ElementClassProperty<element_type>::parent_element_type;
    UInt nb_nodes_parent_el =
        ElementClass<parent_type>::getNbNodesPerInterpolationElement();
    for (UInt i = 0; i < nb_nodes_parent_el; ++i) {
      Vector<Real> parent_c(parent_coords(i));
      parent_c = element_coords(i);
    }
  }

  /// map the points from the reference domain of the subelement to the physical
  /// domain
  static void mapToPhysicalDomain(const Matrix<Real> & element_coords,
                                  Matrix<Real> & sub_coords,
                                  Matrix<Real> & sub_shapes,
                                  Matrix<Real> & physical_points,
                                  UInt sub_element = 0) {
    /// get the sub_element_type

    getSubElementCoords(element_coords, sub_coords, sub_element);
    /// map the points of the subelements in the physical domain
    switch (sub_element) {
    case 0:
      ElementClass<ElementClassProperty<element_type>::sub_element_type_1>::
          interpolate(sub_coords, sub_shapes, physical_points);
      break;
    case 1:
      ElementClass<ElementClassProperty<element_type>::sub_element_type_2>::
          interpolate(sub_coords, sub_shapes, physical_points);
      break;
    }
  }

  /// map the points from the physical domain to the parent reference domain
  static void mapToParentRefDomain(const Matrix<Real> & element_coords,
                                   Matrix<Real> & parent_coords,
                                   Matrix<Real> & physical_points,
                                   Matrix<Real> & natural_coords) {
    const ElementType parent_type =
        ElementClassProperty<element_type>::parent_element_type;
    getParentCoords(element_coords, parent_coords);

    /// map the points from the physical domain into the parent reference domain
    ElementClass<parent_type>::inverseMap(physical_points, parent_coords,
                                          natural_coords);
  }

  /// map the points from the subelement reference domain to the parent
  /// reference domain
  static void mapFromSubRefToParentRef(const Matrix<Real> & element_coords,
                                       Matrix<Real> & parent_coords,
                                       Matrix<Real> & sub_coords,
                                       Matrix<Real> & sub_shapes,
                                       Matrix<Real> & physical_points,
                                       Matrix<Real> & natural_points,
                                       UInt /*nb_points*/, UInt sub_element) {
    mapToPhysicalDomain(element_coords, sub_coords, sub_shapes, physical_points,
                        sub_element);
    mapToParentRefDomain(element_coords, parent_coords, physical_points,
                         natural_points);
  }

  static void mapFromSubRefToParentRef(const Matrix<Real> & element_coords,
                                       Matrix<Real> & sub_coords,
                                       Matrix<Real> & parent_coords,
                                       Matrix<Real> & sub_shapes,
                                       Matrix<Real> & physical_points,
                                       Matrix<Real> & parent_el_natural_coords,
                                       UInt sub_element) {
    mapToPhysicalDomain(element_coords, sub_coords, sub_shapes, physical_points,
                        sub_element);
    mapToParentRefDomain(element_coords, parent_coords, physical_points,
                         parent_el_natural_coords);
  }

  /// compute the normal of a surface defined by the function f
  static inline void
  computeNormalsOnNaturalCoordinates(const Matrix<Real> & /*coord*/,
                                     Matrix<Real> & /*f*/,
                                     Matrix<Real> & /*normals*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /// determine orientation of the element with respect to the interface
  static inline UInt getOrientation(const Vector<bool> & is_inside);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  static AKANTU_GET_MACRO_NOT_CONST(Kind, _ek_igfem, ElementKind);
  static ElementType getP1ElementType() { AKANTU_TO_IMPLEMENT(); };
  static AKANTU_GET_MACRO_NOT_CONST(
      SpatialDimension, ElementClassProperty<element_type>::spatial_dimension,
      UInt);
  static ElementType & getFacetType(UInt /*t*/ = 0) { AKANTU_TO_IMPLEMENT(); }
  static ElementType * getFacetTypeInternal() { AKANTU_TO_IMPLEMENT(); }
private:
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "geometrical_element_igfem.hh"
/* -------------------------------------------------------------------------- */
#include "element_class_igfem_segment_3_inline_impl.hh"
#include "element_class_igfem_triangle_4_inline_impl.hh"
#include "element_class_igfem_triangle_5_inline_impl.hh"
/* -------------------------------------------------------------------------- */

#endif /* AKANTU_ELEMENT_CLASS_IGFEM_HH_ */
