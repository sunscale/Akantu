/**
 * @file   shape_igfem_inline_impl.cc
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Sun May 25 21:53:21 2014
 *
 * @brief  ShapeIGFEM inline implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
ShapeLagrange<_ek_igfem>::ShapeLagrange(const Mesh & mesh,
					const ID & id,
					const MemoryID & memory_id) :
  ShapeFunctions(mesh, id, memory_id),
  shapes("shapes_generic", id),
  shapes_derivatives("shapes_derivatives_generic", id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline const Array<Real> & ShapeLagrange<_ek_igfem>::getShapes(const ElementType & el_type,
							       const GhostType & ghost_type) const {
  return shapes(FEEngine::getInterpolationType(el_type), ghost_type);
}

/* -------------------------------------------------------------------------- */
inline const Array<Real> & ShapeLagrange<_ek_igfem>::getShapesDerivatives(const ElementType & el_type,
									  const GhostType & ghost_type) const {
  return shapes_derivatives(FEEngine::getInterpolationType(el_type), ghost_type);
}

/* -------------------------------------------------------------------------- */
template <ElementType type, bool is_sub>
void ShapeLagrange<_ek_igfem>::setControlPointsByType(const Matrix<Real> & points,
						      const GhostType & ghost_type) {
  if (!is_sub)
    ShapeLagrange<_ek_igfem>::control_points(type, ghost_type).shallowCopy(points);
  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
}

/* -------------------------------------------------------------------------- */
#define INIT_SHAPE_FUNCTIONS(type)					\
  setControlPointsByType<type,ElementClassProperty<type>::is_subelement>(control_points, ghost_type); \
  precomputeShapesOnControlPoints<type,ElementClassProperty<type>::is_subelement>(nodes, ghost_type); \
  if (ElementClass<type>::getNaturalSpaceDimension() ==			\
      mesh.getSpatialDimension() )		\
    precomputeShapeDerivativesOnControlPoints<type,ElementClassProperty<type>::is_subelement>(nodes, ghost_type);

inline void
ShapeLagrange<_ek_igfem>::initShapeFunctions(const Array<Real> & nodes,
					const Matrix<Real> & control_points,
					const ElementType & type,
					const GhostType & ghost_type) {
  AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(INIT_SHAPE_FUNCTIONS);
}

#undef INIT_SHAPE_FUNCTIONS

/* -------------------------------------------------------------------------- */

template <ElementType type, bool is_sub>
inline void ShapeLagrange<_ek_igfem>::
computeShapeDerivativesOnCPointsByElement(const Matrix<Real> & node_coords,
					  const Matrix<Real> & natural_coords,
					  Tensor3<Real> & shapesd) {

  if(!is_sub) {
    // 
    // compute dnds
    Tensor3<Real> dnds(node_coords.rows(), node_coords.cols(), natural_coords.cols());
    ElementClass<type>::computeDNDS(natural_coords, dnds);
    // compute dxds
    Tensor3<Real> J(node_coords.rows(), natural_coords.rows(), natural_coords.cols());
    ElementClass<type>::computeJMat(dnds, node_coords, J);

    // compute shape derivatives
    ElementClass<type>::computeShapeDerivatives(J, dnds, shapesd); 

  }

  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_igfem>::inverseMap(const Vector<Real> & real_coords,
					  UInt elem,
					  Vector<Real> & natural_coords,
					  const GhostType & ghost_type) const{
  inverseMap<type, ElementClassProperty<type>::is_subelement>(real_coords, elem, natural_coords, ghost_type);

}

/* -------------------------------------------------------------------------- */
template <ElementType type, bool is_sub>
void ShapeLagrange<_ek_igfem>::inverseMap(const Vector<Real> & real_coords,
					  UInt elem,
					  Vector<Real> & natural_coords,
					  const GhostType & ghost_type) const{

  if (!is_sub) {
    UInt spatial_dimension = mesh.getSpatialDimension();
    UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerInterpolationElement();

    UInt * elem_val = mesh.getConnectivity(type, ghost_type).storage();
    Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

    mesh.extractNodalValuesFromElement(mesh.getNodes(),
				       nodes_coord.storage(),
				       elem_val + elem*nb_nodes_per_element,
				       nb_nodes_per_element,
				       spatial_dimension);

    ElementClass<type>::inverseMap(real_coords,
				   nodes_coord,
				   natural_coords);
  }

  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
bool ShapeLagrange<_ek_igfem>::contains(const Vector<Real> & real_coords,
					UInt elem,
					const GhostType & ghost_type) const{

  return contains<type,ElementClassProperty<type>::is_subelement>(real_coords, elem, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <ElementType type,bool is_sub>
bool ShapeLagrange<_ek_igfem>::contains(const Vector<Real> & real_coords,
					UInt elem,
					const GhostType & ghost_type) const{

  if (!is_sub) {
    UInt spatial_dimension = mesh.getSpatialDimension();
    Vector<Real> natural_coords(spatial_dimension);

    inverseMap<type>(real_coords, elem, natural_coords, ghost_type);
    return ElementClass<type>::contains(natural_coords);
  }
  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_igfem>::computeShapes(const Vector<Real> & real_coords,
					     UInt elem,
					     Vector<Real> & shapes,
					     const GhostType & ghost_type) const{

  computeShapes<type,ElementClassProperty<type>::is_subelement>(real_coords, elem, shapes, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <ElementType type, bool is_sub>
void ShapeLagrange<_ek_igfem>::computeShapes(const Vector<Real> & real_coords,
					     UInt elem,
					     Vector<Real> & shapes,
					     const GhostType & ghost_type) const{

  if (!is_sub) {
    UInt spatial_dimension = mesh.getSpatialDimension();
    Vector<Real> natural_coords(spatial_dimension);

    inverseMap<type>(real_coords, elem, natural_coords, ghost_type);
    ElementClass<type>::computeShapes(natural_coords, shapes);
  }
  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, bool is_sub>
void ShapeLagrange<_ek_igfem>::precomputeShapesOnControlPoints(__attribute__((unused)) const Array<Real> & nodes,
							       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if (!is_sub) {
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

  }

  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type, bool is_sub>
void ShapeLagrange<_ek_igfem>::precomputeShapeDerivativesOnControlPoints(const Array<Real> & nodes, GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  if (!is_sub) {

    InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

    //  Real * coord = mesh.getNodes().storage();
    UInt spatial_dimension = mesh.getSpatialDimension();

    UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerInterpolationElement();
    UInt size_of_shapesd      = ElementClass<type>::getShapeDerivativesSize();
    Matrix<Real> & natural_coords = control_points(type, ghost_type);
    UInt nb_points = natural_coords.cols();

    UInt nb_element = mesh.getConnectivity(type, ghost_type).getSize();

    Array<Real> & shapes_derivatives_tmp = shapes_derivatives.alloc(nb_element*nb_points,
								    size_of_shapesd,
								    itp_type,
								    ghost_type);

    Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
    FEEngine::extractNodalToElementField(mesh, nodes, x_el,
				    type, ghost_type);

    Real * shapesd_val = shapes_derivatives_tmp.storage();
    Array<Real>::matrix_iterator x_it = x_el.begin(spatial_dimension,
						   nb_nodes_per_element);

    for (UInt elem = 0; elem < nb_element; ++elem, ++x_it) {

      Matrix<Real> & X = *x_it;
      Tensor3<Real> B(shapesd_val,
		      spatial_dimension, nb_nodes_per_element, nb_points);


      computeShapeDerivativesOnCPointsByElement<type, ElementClassProperty<type>::is_subelement>(X,
												 natural_coords,
												 B);

      shapesd_val += size_of_shapesd*nb_points;
    }
  }

  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_igfem>::interpolateOnControlPoints(const Array<Real> &in_u,
							  Array<Real> &out_uq,
							  UInt nb_degree_of_freedom,
							  GhostType ghost_type,
							  const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();
  interpolateOnControlPoints<type, ElementClassProperty<type>::is_subelement>(in_u, out_uq, 
									      nb_degree_of_freedom, 
									      ghost_type, filter_elements);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type, bool is_sub>
void ShapeLagrange<_ek_igfem>::interpolateOnControlPoints(const Array<Real> &in_u,
							  Array<Real> &out_uq,
							  UInt nb_degree_of_freedom,
							  GhostType ghost_type,
							  const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  if(!is_sub) {
    InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
    AKANTU_DEBUG_ASSERT(shapes.exists(itp_type, ghost_type),
			"No shapes for the type "
			<< shapes.printType(itp_type, ghost_type));

    UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerInterpolationElement();

    Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
    FEEngine::extractNodalToElementField(mesh, in_u, u_el, type, ghost_type, filter_elements);

    this->interpolateElementalFieldOnControlPoints<type>(u_el, out_uq, ghost_type,
							 shapes(itp_type, ghost_type),
							 filter_elements);
  }
  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_igfem>::gradientOnControlPoints(const Array<Real> &in_u,
						       Array<Real> &out_nablauq,
						       UInt nb_degree_of_freedom,
						       GhostType ghost_type,
						       const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  gradientOnControlPoints<type,ElementClassProperty<type>::is_subelement>(in_u, out_nablauq, nb_degree_of_freedom, 
									  ghost_type, filter_elements);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type, bool is_sub>
void ShapeLagrange<_ek_igfem>::gradientOnControlPoints(const Array<Real> &in_u,
						       Array<Real> &out_nablauq,
						       UInt nb_degree_of_freedom,
						       GhostType ghost_type,
						       const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  if (!is_sub) {
    InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
    AKANTU_DEBUG_ASSERT(shapes_derivatives.exists(itp_type, ghost_type),
			"No shapes derivatives for the type "
			<< shapes_derivatives.printType(itp_type, ghost_type));

    UInt nb_nodes_per_element  = ElementClass<type>::getNbNodesPerInterpolationElement();

    Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
    FEEngine::extractNodalToElementField(mesh, in_u, u_el, type, ghost_type, filter_elements);

    this->gradientElementalFieldOnControlPoints<type>(u_el, out_nablauq, ghost_type,
						      shapes_derivatives(itp_type, ghost_type),
						      filter_elements);
  }

  else {

    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_igfem>::fieldTimesShapes(const Array<Real> & field,
						Array<Real> & field_times_shapes,
						GhostType ghost_type) const {

  fieldTimesShapes<type, ElementClassProperty<type>::is_subelement>(field, field_times_shapes, ghost_type);

}

/* -------------------------------------------------------------------------- */
template <ElementType type, bool is_sub>
void ShapeLagrange<_ek_igfem>::fieldTimesShapes(const Array<Real> & field,
						Array<Real> & field_times_shapes,
						GhostType ghost_type) const {

  if (!is_sub) {
    field_times_shapes.resize(field.getSize());

    UInt size_of_shapes = ElementClass<type>::getShapeSize();
    InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
    UInt nb_degree_of_freedom = field.getNbComponent();

    const Array<Real> & shape = shapes(itp_type, ghost_type);

    Array<Real>::const_matrix_iterator field_it  = field.begin(nb_degree_of_freedom, 1);
    Array<Real>::const_matrix_iterator shapes_it = shape.begin(1, size_of_shapes);

    Array<Real>::matrix_iterator it  = field_times_shapes.begin(nb_degree_of_freedom, size_of_shapes);
    Array<Real>::matrix_iterator end = field_times_shapes.end  (nb_degree_of_freedom, size_of_shapes);

    for (; it != end; ++it, ++field_it, ++shapes_it) {
      it->mul<false, false>(*field_it, *shapes_it);
    }
  }
  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
}

/* -------------------------------------------------------------------------- */
void ShapeLagrange<_ek_igfem>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Shapes Lagrange [" << std::endl;
  ShapeFunctions::printself(stream, indent + 1);
  shapes.printself(stream, indent + 1);
  shapes_derivatives.printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
template <bool is_sub>
void ShapeLagrange<_ek_igfem>::printself(std::ostream & stream, int indent) const {

  if (!is_sub) {
    std::string space;
    for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

    stream << space << "Shapes Lagrange [" << std::endl;
    ShapeFunctions::printself(stream, indent + 1);
    shapes.printself(stream, indent + 1);
    shapes_derivatives.printself(stream, indent + 1);
    stream << space << "]" << std::endl;
  }
  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
}

/* -------------------------------------------------------------------------- */


// __END_AKANTU__

// #include "fe_engine.hh"

// __BEGIN_AKANTU__

// /* -------------------------------------------------------------------------- */
// template <ElementKind kind>
// template <ElementType type>
// void ShapeLagrange<kind>::precomputeShapesOnControlPoints(__attribute__((unused)) const Array<Real> & nodes,
// 							  GhostType ghost_type) {
//   AKANTU_DEBUG_IN();

//   InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

//   /// get natural coordinates of subelement
//   Matrix<Real> & natural_coords = control_points(type, ghost_type);
//   UInt nb_points = natural_coords.cols();


//   UInt size_of_shapes = ElementClass<type>::getShapeSize();

//   UInt nb_element = mesh.getConnectivity(type, ghost_type).getSize();;

//   Array<Real> & shapes_tmp = shapes.alloc(nb_element*nb_points,
// 					   size_of_shapes,
// 					   itp_type,
// 					   ghost_type);

//   Array<Real>::matrix_iterator shapes_it =
//     shapes_tmp.begin_reinterpret(ElementClass<type>::getNbNodesPerInterpolationElement(), nb_points,
// 				 nb_element);

//   Array<Real> & parent_control_points_tmp = parent_control_points.alloc(nb_element,
// 									nb_points,
// 									itp_type,
// 									ghost_type);

//   Array<Real>::matrix_iterator parent_control_points_it =
//     parent_control_points_tmp.begin(nb_points,spatial_dimension);

//   Matrix<Real> physical_control_points(nb_points,nb_points);
//   physical_control_points.clear();

//   Matrix<Real> N_tmp = ();
//   for (UInt elem = 0; elem < nb_element; ++elem, ++shapes_it; ++parent_control_points_it) {
//     Matrix<Real> & N = *shapes_it;
//     ElementClass<type>::computeShapes(natural_coords,
// 				      N_tmp);

//     if (!ElementClassProperty<type>::has_parent)
//       N = N_tmp;
//     else {
//       /// compute control points in physical domain
//       physical_control_points.mul<false, false>(natural_coords, control_points);

//       Array<Real>::vector_iterator physical_control_points_it = physical_control_points.begin(spatial_dimension);

//       /// compute control points in reference domain of the parent 
//       for (UInt i = 0; i < spatial_dimension; ++i) {
// 	inverseMap(physical_control_points, elem, *physical_control_points_it, ghost_type);
//       }

//       /// compute parent shape functions
//       ElementClass<ElementClassProperty<type>::parent_el_type>::computeShapes(natural_coords_parents)

//     }
//   }
//   AKANTU_DEBUG_OUT();
// }

// /* -------------------------------------------------------------------------- */
// template <ElementKind kind>
// template <ElementType type>
// void ShapeLagrange<kind>::inverseMap(const Vector<Real> & real_coords,
// 				     UInt elem,
// 				     Vector<Real> & natural_coords,
// 				     const GhostType & ghost_type) const{

//   UInt spatial_dimension = mesh.getSpatialDimension();
//   UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerInterpolationElement();

//   UInt * elem_val = mesh.getConnectivity(type, ghost_type).storage();
//   Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

//   mesh.extractNodalValuesFromElement(mesh.getNodes(),
// 				     nodes_coord.storage(),
// 				     elem_val + elem*nb_nodes_per_element,
// 				     nb_nodes_per_element,
// 				     spatial_dimension);

//   ElementClass<type>::inverseMap(real_coords,
// 				 nodes_coord,
// 				 natural_coords);
// }

// /* -------------------------------------------------------------------------- */
