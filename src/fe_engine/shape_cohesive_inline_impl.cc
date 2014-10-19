/**
 * @file   shape_cohesive_inline_impl.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 23 2012
 * @date last modification: Fri Jun 13 2014
 *
 * @brief  ShapeCohesive inline implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
inline ShapeLagrange<_ek_cohesive>::ShapeLagrange(const Mesh & mesh,
						  const ID & id,
						  const MemoryID & memory_id) :
  ShapeFunctions(mesh, id, memory_id),
  shapes("shapes_cohesive", id),
  shapes_derivatives("shapes_derivatives_cohesive", id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

#define INIT_SHAPE_FUNCTIONS(type)					\
  setControlPointsByType<type>(control_points, ghost_type);		\
  precomputeShapesOnControlPoints<type>(nodes, ghost_type);		\
  precomputeShapeDerivativesOnControlPoints<type>(nodes, ghost_type);

/* -------------------------------------------------------------------------- */
inline void ShapeLagrange<_ek_cohesive>::initShapeFunctions(const Array<Real> & nodes,
						const Matrix<Real> & control_points,
						const ElementType & type,
						const GhostType & ghost_type) {
  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(INIT_SHAPE_FUNCTIONS);
}

/* -------------------------------------------------------------------------- */
inline const Array<Real> & ShapeLagrange<_ek_cohesive>::getShapes(const ElementType & el_type,
								   const GhostType & ghost_type) const {
  return shapes(FEEngine::getInterpolationType(el_type), ghost_type);
}

/* -------------------------------------------------------------------------- */
inline const Array<Real> & ShapeLagrange<_ek_cohesive>::getShapesDerivatives(const ElementType & el_type,
									      const GhostType & ghost_type) const {
  return shapes_derivatives(FEEngine::getInterpolationType(el_type), ghost_type);
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_cohesive>::precomputeShapesOnControlPoints(__attribute__((unused)) const Array<Real> & nodes,
								  GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

  Matrix<Real> & natural_coords = control_points(type, ghost_type);
  UInt nb_points = natural_coords.cols();

  UInt size_of_shapes = ElementClass<type>::getShapeSize();

  UInt nb_element = mesh.getConnectivity(type, ghost_type).getSize();;

  Array<Real> & shapes_tmp = shapes.alloc(nb_element*nb_points,
					   size_of_shapes,
					   itp_type,
					   ghost_type);

  Array<Real>::matrix_iterator shapes_it =
    shapes_tmp.begin_reinterpret(ElementClass<type>::getNbNodesPerInterpolationElement(), nb_points,
				 nb_element);

  for (UInt elem = 0; elem < nb_element; ++elem, ++shapes_it) {
    Matrix<Real> & N = *shapes_it;
    ElementClass<type>::computeShapes(natural_coords,
				      N);
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_cohesive>::precomputeShapeDerivativesOnControlPoints(__attribute__((unused)) const Array<Real> & nodes,
									    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt size_of_shapesd      = ElementClass<type>::getShapeDerivativesSize();
  UInt spatial_dimension    = ElementClass<type>::getNaturalSpaceDimension();
  UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerInterpolationElement();

  Matrix<Real> natural_coords = this->control_points(type, ghost_type);
  UInt nb_points = natural_coords.cols();

  // UInt * elem_val = this->mesh->getConnectivity(type, ghost_type).storage();;
  UInt nb_element = this->mesh.getConnectivity(type, ghost_type).getSize();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

  Array<Real> & shapes_derivatives_tmp =
    this->shapes_derivatives.alloc(nb_element*nb_points,
				   size_of_shapesd,
				   itp_type,
				   ghost_type);

  Real * shapesd_val = shapes_derivatives_tmp.storage();
  for (UInt elem = 0; elem < nb_element; ++elem) {
    Tensor3<Real> B(shapesd_val,
			   spatial_dimension, nb_nodes_per_element, nb_points);
    ElementClass<type>::computeDNDS(natural_coords, B);

    shapesd_val += size_of_shapesd*nb_points;
  }

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::extractNodalToElementField(const Array<Real> & nodal_f,
							     Array<Real> & elemental_f,
							     const GhostType & ghost_type,
							     const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_itp_element = ElementClass<type>::getNbNodesPerInterpolationElement();
  UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();
  UInt nb_degree_of_freedom = nodal_f.getNbComponent();
  UInt nb_element = this->mesh.getNbElement(type, ghost_type);
  UInt * conn_val = this->mesh.getConnectivity(type, ghost_type).storage();

  if(filter_elements != empty_filter) {
    nb_element      = filter_elements.getSize();
  }

  elemental_f.resize(nb_element);

  Array<Real>::matrix_iterator u_it = elemental_f.begin(nb_degree_of_freedom,
								  nb_nodes_per_itp_element);

  UInt * el_conn;
  ReduceFunction reduce_function;

  for (UInt el = 0; el < nb_element; ++el, ++u_it) {
    Matrix<Real> & u = *u_it;
    if(filter_elements != empty_filter) el_conn = conn_val + filter_elements(el) * nb_nodes_per_element;
    else el_conn = conn_val + el * nb_nodes_per_element;

    // compute the average/difference of the nodal field loaded from cohesive element
    for (UInt n = 0; n < nb_nodes_per_itp_element; ++n) {
      UInt node_plus = *(el_conn + n);
      UInt node_minus = *(el_conn + n + nb_nodes_per_itp_element);
      for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
	Real u_plus  = nodal_f(node_plus, d);
	Real u_minus = nodal_f(node_minus, d);
	u(d, n)	= reduce_function(u_plus, u_minus);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::interpolateOnControlPoints(const Array<Real> &in_u,
							     Array<Real> &out_uq,
							     UInt nb_degree_of_freedom,
							     GhostType ghost_type,
							     const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

  AKANTU_DEBUG_ASSERT(this->shapes.exists(itp_type, ghost_type),
		      "No shapes for the type "
		      << this->shapes.printType(itp_type, ghost_type));

  UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerInterpolationElement();
  Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  this->extractNodalToElementField<type, ReduceFunction>(in_u, u_el, ghost_type, filter_elements);

  this->template interpolateElementalFieldOnControlPoints<type>(u_el, out_uq, ghost_type,
								shapes(itp_type, ghost_type),
								filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::variationOnControlPoints(const Array<Real> &in_u,
							   Array<Real> &nablauq,
							   UInt nb_degree_of_freedom,
							   GhostType ghost_type,
							   const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

  AKANTU_DEBUG_ASSERT(this->shapes_derivatives.exists(itp_type, ghost_type),
		      "No shapes for the type "
		      << this->shapes_derivatives.printType(itp_type, ghost_type));

  UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerInterpolationElement();
  Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  this->extractNodalToElementField<type, ReduceFunction>(in_u, u_el, ghost_type, filter_elements);

  this->template gradientElementalFieldOnControlPoints<type>(u_el, nablauq, ghost_type,
							     shapes_derivatives(itp_type, ghost_type),
							     filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <ElementType type, class ReduceFunction>
void ShapeLagrange<_ek_cohesive>::computeNormalsOnControlPoints(const Array<Real> &u,
								Array<Real> &normals_u,
								GhostType ghost_type,
								const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_element = this->mesh.getNbElement(type, ghost_type);
  UInt nb_points  = this->control_points(type, ghost_type).cols();
  UInt spatial_dimension = this->mesh.getSpatialDimension();

  if(filter_elements != empty_filter)
    nb_element  = filter_elements.getSize();

  normals_u.resize(nb_points * nb_element);

  Array<Real> tangents_u(nb_element * nb_points,
			 (spatial_dimension *  (spatial_dimension -1)));

  this->template variationOnControlPoints<type, ReduceFunction>(u,
								tangents_u,
								spatial_dimension,
								ghost_type,
								filter_elements);

  Array<Real>::vector_iterator normal     = normals_u.begin(spatial_dimension);
  Array<Real>::vector_iterator normal_end = normals_u.end(spatial_dimension);

  Real * tangent = tangents_u.storage();

  if(spatial_dimension == 3)
    for (; normal != normal_end; ++normal) {
      Math::vectorProduct3(tangent, tangent+spatial_dimension, normal->storage());

      (*normal) /= normal->norm();

      tangent += spatial_dimension * 2;
    }
  else if (spatial_dimension == 2)
    for (; normal != normal_end; ++normal) {
      Vector<Real> a1(tangent, spatial_dimension);

      (*normal)(0) = -a1(1);
      (*normal)(1) = a1(0);
      (*normal) /= normal->norm();

      tangent += spatial_dimension;
    }

  AKANTU_DEBUG_OUT();
}

