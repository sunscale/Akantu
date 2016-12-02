/**
 * @file   shape_linked_inline_impl.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  ShapeLinked inline implementation
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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

template <ElementKind kind>
inline void ShapeLinked<kind>::initShapeFunctions(
    __attribute__((unused)) const Array<Real> & nodes,
    __attribute__((unused)) const Matrix<Real> & integration_points,
    __attribute__((unused)) const ElementType & type,
    __attribute__((unused)) const GhostType & ghost_type) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

#undef INIT_SHAPE_FUNCTIONS

/* -------------------------------------------------------------------------- */
#define INIT_SHAPE_FUNCTIONS(type)                                             \
  setIntegrationPointsByType<type>(integration_points, ghost_type);            \
  precomputeShapesOnIntegrationPoints<type>(nodes, ghost_type);                \
  precomputeShapeDerivativesOnIntegrationPoints<type>(nodes, ghost_type);

#if defined(AKANTU_STRUCTURAL_MECHANICS)
template <>
inline void ShapeLinked<_ek_structural>::initShapeFunctions(
    __attribute__((unused)) const Array<Real> & nodes,
    __attribute__((unused)) const Matrix<Real> & integration_points,
    __attribute__((unused)) const ElementType & type,
    __attribute__((unused)) const GhostType & ghost_type) {
  AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(INIT_SHAPE_FUNCTIONS);
}
#endif

#undef INIT_SHAPE_FUNCTIONS

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
inline const Array<Real> &
ShapeLinked<kind>::getShapes(const ElementType & type,
                             const GhostType & ghost_type, UInt id) const {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(shapes.exists(type, ghost_type),
                      "No shapes of type " << type << " in " << this->id);
  AKANTU_DEBUG_OUT();
  return *(shapes(type, ghost_type)[id]);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
inline const Array<Real> & ShapeLinked<kind>::getShapesDerivatives(
    const ElementType & type, const GhostType & ghost_type, UInt id) const {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(shapes_derivatives.exists(type, ghost_type),
                      "No shapes_derivatives of type " << type << " in "
                                                       << this->id);
  AKANTU_DEBUG_OUT();
  return *(shapes_derivatives(type, ghost_type)[id]);
}

#if defined(AKANTU_STRUCTURAL_MECHANICS)
/* -------------------------------------------------------------------------- */
template <>
template <ElementType type>
void ShapeLinked<_ek_structural>::precomputeShapesOnIntegrationPoints(
    const Array<Real> & nodes, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  //  Real * coord = mesh.getNodes().storage();
  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  Matrix<Real> & natural_coords = integration_points(type, ghost_type);
  UInt nb_points = integration_points(type, ghost_type).cols();

  UInt size_of_shapes = ElementClass<type>::getShapeSize();

  std::string ghost = "";
  if (ghost_type == _ghost) {
    ghost = "ghost_";
  }

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  UInt nb_shape_functions =
      ElementClass<type, _ek_structural>::getNbShapeFunctions();

  Array<Real> ** shapes_tmp = new Array<Real> * [nb_shape_functions];

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type);

  for (UInt s = 0; s < nb_shape_functions; ++s) {
    std::stringstream sstr_shapes;
    sstr_shapes << id << ":" << ghost << "shapes:" << type << ":" << s;
    shapes_tmp[s] = &(alloc<Real>(sstr_shapes.str(), nb_element * nb_points,
                                  size_of_shapes));
    Array<Real>::matrix_iterator x_it =
        x_el.begin(spatial_dimension, nb_nodes_per_element);
    Array<Real>::matrix_iterator shapes_it =
        shapes_tmp[s]->begin_reinterpret(size_of_shapes, nb_points, nb_element);

    for (UInt elem = 0; elem < nb_element; ++elem, ++shapes_it, ++x_it) {
      Matrix<Real> & X = *x_it;
      Matrix<Real> & N = *shapes_it;
      ElementClass<type>::computeShapes(natural_coords, N, X, s);
    }
  }

  shapes(type, ghost_type) = shapes_tmp;

  AKANTU_DEBUG_OUT();
}
#endif

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLinked<kind>::precomputeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  // Real * coord = mesh.getNodes().storage();
  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt natural_spatial_dimension =
      ElementClass<type>::getNaturalSpaceDimension();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt size_of_shapesd = ElementClass<type>::getShapeDerivativesSize();
  Matrix<Real> & natural_coords = integration_points(type, ghost_type);
  UInt nb_points = natural_coords.cols();

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  std::string ghost = "";

  if (ghost_type == _ghost) {
    ghost = "ghost_";
  }

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type);

  UInt nb_shape_functions = ElementClass<type>::getNbShapeDerivatives();

  Array<Real> ** shapes_derivatives_tmp =
      new Array<Real> * [nb_shape_functions];
  for (UInt s = 0; s < nb_shape_functions; ++s) {
    std::stringstream sstr_shapesd;
    sstr_shapesd << id << ":" << ghost << "shapes_derivatives:" << type << ":"
                 << s;
    shapes_derivatives_tmp[s] = &(alloc<Real>(
        sstr_shapesd.str(), nb_element * nb_points, size_of_shapesd));
    Real * shapesd_val = shapes_derivatives_tmp[s]->storage();

    Array<Real>::matrix_iterator x_it =
        x_el.begin(spatial_dimension, nb_nodes_per_element);

    for (UInt elem = 0; elem < nb_element; ++elem, ++x_it) {
      // compute shape derivatives
      Matrix<Real> & X = *x_it;
      Tensor3<Real> B(shapesd_val, natural_spatial_dimension,
                      nb_nodes_per_element, nb_points);
      ElementClass<type>::computeShapeDerivatives(natural_coords, B, X, s);

      shapesd_val += size_of_shapesd * nb_points;
    }
  }

  shapes_derivatives(type, ghost_type) = shapes_derivatives_tmp;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLinked<kind>::extractNodalToElementField(
    const Array<Real> & nodal_f, Array<Real> & elemental_f,
    UInt num_degre_of_freedom_to_extract, const GhostType & ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_degree_of_freedom = nodal_f.getNbComponent();
  UInt nb_element = mesh.getNbElement(type, ghost_type);
  UInt * conn_val = mesh.getConnectivity(type, ghost_type).storage();

  if (filter_elements != empty_filter) {
    nb_element = filter_elements.getSize();
  }

  elemental_f.resize(nb_element);

  Real * nodal_f_val = nodal_f.storage();
  Real * f_val = elemental_f.storage();

  UInt * el_conn;
  for (UInt el = 0; el < nb_element; ++el) {
    if (filter_elements != empty_filter)
      el_conn = conn_val + filter_elements(el) * nb_nodes_per_element;
    else
      el_conn = conn_val + el * nb_nodes_per_element;

    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt node = *(el_conn + n);
      *f_val = nodal_f_val[node * nb_degree_of_freedom +
                           num_degre_of_freedom_to_extract];
      f_val += 1;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLinked<kind>::interpolateOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_uq,
    __attribute__((unused)) UInt nb_degree_of_freedom,
    const GhostType & ghost_type, const Array<UInt> & filter_elements,
    bool accumulate, UInt id_shape, UInt num_degre_of_freedom_to_interpolate,
    __attribute__((unused)) UInt num_degre_of_freedom_interpolated) const {
  AKANTU_DEBUG_IN();

  Array<Real> * shapes_loc = shapes(type, ghost_type)[id_shape];

  AKANTU_DEBUG_ASSERT(shapes_loc != NULL, "No shapes for the type " << type);

  UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();
  Array<Real> u_el(0, nb_nodes_per_element);
  extractNodalToElementField<type>(in_u, u_el,
                                   num_degre_of_freedom_to_interpolate,
                                   ghost_type, filter_elements);

  if (!accumulate)
    out_uq.clear();

  UInt nb_points = integration_points(type, ghost_type).cols() * u_el.getSize();
  Array<Real> uq(nb_points, 1, 0.);

  this->template interpolateElementalFieldOnIntegrationPoints<type>(
      u_el, uq, ghost_type, *shapes_loc, filter_elements);

  for (UInt q = 0; q < nb_points; ++q) {
    out_uq(q, num_degre_of_freedom_to_interpolate) += uq(q);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLinked<kind>::gradientOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_nablauq,
    UInt nb_degree_of_freedom, const GhostType & ghost_type,
    const Array<UInt> & filter_elements, bool accumulate, UInt id_shape,
    UInt num_degre_of_freedom_to_interpolate,
    __attribute__((unused)) UInt num_degre_of_freedom_interpolated) const {
  AKANTU_DEBUG_IN();

  Array<Real> * shapesd_loc = shapes_derivatives(type, ghost_type)[id_shape];

  AKANTU_DEBUG_ASSERT(shapesd_loc != NULL, "No shapes for the type " << type);

  UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();
  Array<Real> u_el(0, nb_nodes_per_element);
  extractNodalToElementField<type>(in_u, u_el,
                                   num_degre_of_freedom_to_interpolate,
                                   ghost_type, filter_elements);

  UInt nb_points = integration_points(type, ghost_type).cols() * u_el.getSize();
  UInt element_dimension = ElementClass<type>::getSpatialDimension();

  Array<Real> nablauq(nb_points, element_dimension, 0.);

  if (!accumulate)
    out_nablauq.clear();
  this->template gradientElementalFieldOnIntegrationPoints<type>(
      u_el, nablauq, ghost_type, *shapesd_loc, filter_elements);

  Array<Real>::matrix_iterator nabla_u_it = nablauq.begin(1, element_dimension);
  Array<Real>::matrix_iterator out_nabla_u_it =
      out_nablauq.begin(nb_degree_of_freedom, element_dimension);
  for (UInt q = 0; q < nb_points; ++q, ++nabla_u_it, ++out_nabla_u_it) {
    for (UInt s = 0; s < element_dimension; ++s) {
      (*out_nabla_u_it)(num_degre_of_freedom_to_interpolate, s) +=
          (*nabla_u_it)(0, s);
    }
  }

  AKANTU_DEBUG_OUT();
}
