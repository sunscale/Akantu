/**
 * @file   shape_igfem_inline_impl.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  ShapeIGFEM inline implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

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
#define INIT_SHAPE_FUNCTIONS(type)					\
  setControlPointsByType<type>(control_points, ghost_type);		\
  setControlPointsByType<ElementClassProperty<type>::sub_element_type_1>(control_points_1, ghost_type); \
  setControlPointsByType<ElementClassProperty<type>::sub_element_type_2>(control_points_2, ghost_type); \
  precomputeShapesOnControlPoints<type>(nodes, ghost_type);		\
  if (ElementClass<type>::getNaturalSpaceDimension() ==			\
      mesh.getSpatialDimension())					\
    precomputeShapeDerivativesOnControlPoints<type>(nodes, ghost_type);


inline void ShapeLagrange<_ek_igfem>::initShapeFunctions(const Array<Real> & nodes,
							 const Matrix<Real> & control_points, 
							 const Matrix<Real> & control_points_1,
							 const Matrix<Real> & control_points_2,
							 const ElementType & type,
							 const GhostType & ghost_type) {
  
  AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(INIT_SHAPE_FUNCTIONS);
}
#undef INIT_SHAPE_FUNCTIONS

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void ShapeLagrange<_ek_igfem>::
computeShapeDerivativesOnCPointsByElement(const Matrix<Real> & node_coords,
					  const Matrix<Real> & natural_coords,
					  Tensor3<Real> & shapesd) const {
  AKANTU_DEBUG_IN();

  // compute dnds
  Tensor3<Real> dnds(node_coords.rows(), node_coords.cols(), natural_coords.cols());
  ElementClass<type>::computeDNDS(natural_coords, dnds);
  // compute dxds
  Tensor3<Real> J(node_coords.rows(), natural_coords.rows(), natural_coords.cols());
  ElementClass<type>::computeJMat(dnds, node_coords, J);

  // compute shape derivatives
  ElementClass<type>::computeShapeDerivatives(J, dnds, shapesd);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_igfem>::inverseMap(const Vector<Real> & real_coords,
				     UInt elem,
				     Vector<Real> & natural_coords,
				     const GhostType & ghost_type) const{

  AKANTU_DEBUG_IN();

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

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
bool ShapeLagrange<_ek_igfem>::contains(const Vector<Real> & real_coords,
			     UInt elem,
			     const GhostType & ghost_type) const{

  UInt spatial_dimension = mesh.getSpatialDimension();
  Vector<Real> natural_coords(spatial_dimension);

  inverseMap<type>(real_coords, elem, natural_coords, ghost_type);
  return ElementClass<type>::contains(natural_coords);
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_igfem>::computeShapes(const Vector<Real> & real_coords,
				  UInt elem,
				  Vector<Real> & shapes,
				  const GhostType & ghost_type) const {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  Vector<Real> natural_coords(spatial_dimension);

  inverseMap<type>(real_coords, elem, natural_coords, ghost_type);
  ElementClass<type>::computeShapes(natural_coords, shapes);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_igfem>::computeShapeDerivatives(const Matrix<Real> & real_coords,
						    UInt elem,
						    Tensor3<Real> & shapesd,
						    const GhostType & ghost_type) const {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_points = real_coords.cols();
  UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerInterpolationElement();

  AKANTU_DEBUG_ASSERT(mesh.getSpatialDimension() == shapesd.size(0) && nb_nodes_per_element == shapesd.size(1),
		      "Shape size doesn't match");
  AKANTU_DEBUG_ASSERT(nb_points == shapesd.size(2),
		      "Number of points doesn't match shapes size");

  Matrix<Real> natural_coords(spatial_dimension, nb_points);

  // Creates the matrix of natural coordinates
  for (UInt i = 0 ; i < nb_points ; i++) {
    Vector<Real> real_point = real_coords(i);
    Vector<Real> natural_point = natural_coords(i);

    inverseMap<type>(real_point, elem, natural_point, ghost_type);
  }


  UInt * elem_val = mesh.getConnectivity(type, ghost_type).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(mesh.getNodes(),
				     nodes_coord.storage(),
				     elem_val + elem*nb_nodes_per_element,
				     nb_nodes_per_element,
				     spatial_dimension);

  computeShapeDerivativesOnCPointsByElement<type>(nodes_coord, natural_coords, shapesd);

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_igfem>::precomputeShapesOnControlPoints(__attribute__((unused)) const Array<Real> & nodes,
							  GhostType ghost_type) {
  AKANTU_DEBUG_IN();
 
  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

  /// typedef for the two subelement_types and the parent element type
  const ElementType sub_type_1 = ElementClassProperty<type>::sub_element_type_1;
  const ElementType sub_type_2 = ElementClassProperty<type>::sub_element_type_2;
  const ElementType parent_type = ElementClassProperty<type>::parent_element_type;

  /// get the spatial dimension for the given element type
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  /// get the integration points for the subelements 
  Matrix<Real> & natural_coords_sub_1 = control_points(sub_type_1, ghost_type);
  Matrix<Real> & natural_coords_sub_2 = control_points(sub_type_2, ghost_type);

  /// store the number of quadrature points on each subelement and the toal number 
  UInt nb_points_sub_1 = natural_coords_sub_1.cols();
  UInt nb_points_sub_2 = natural_coords_sub_2.cols();
  UInt nb_total_points = nb_points_sub_1 + nb_points_sub_2;

  // get the integration points for the parent element
  UInt nb_element = mesh.getConnectivity(type, ghost_type).getSize();
  Array<Real> & natural_coords_parent = igfem_control_points.alloc(nb_element*nb_total_points,
								   spatial_dimension,
								   type,
								   ghost_type);
  Array<Real>::matrix_iterator natural_coords_parent_it = natural_coords_parent.begin_reinterpret(spatial_dimension, nb_total_points, nb_element);
  
  /// get the size of the shapes 
  UInt size_of_shapes = ElementClass<type>::getShapeSize();
  UInt size_of_parent_shapes = ElementClass<parent_type>::getShapeSize();
  UInt size_of_sub_1_shapes = ElementClass<sub_type_1>::getShapeSize();
  UInt size_of_sub_2_shapes = ElementClass<sub_type_2>::getShapeSize();

  /// initialize the matrices to store the shape functions of the subelements and the parent
  Matrix<Real> sub_1_shapes(size_of_sub_1_shapes, nb_points_sub_1);
  Matrix<Real> sub_2_shapes(size_of_sub_2_shapes, nb_points_sub_2);
  Matrix<Real> parent_1_shapes(size_of_parent_shapes, nb_points_sub_1);
  Matrix<Real> parent_2_shapes(size_of_parent_shapes, nb_points_sub_2);
  /// compute the shape functions of the subelements
  ElementClass<sub_type_1>::computeShapes(natural_coords_sub_1, sub_1_shapes);
  ElementClass<sub_type_2>::computeShapes(natural_coords_sub_2, sub_2_shapes);

  /// get the nodal coordinates per element
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el,
				  type, ghost_type);
  Array<Real>::matrix_iterator x_it = x_el.begin(spatial_dimension,
						 nb_nodes_per_element);

  /// allocate the shapes for the given element type
  Array<Real> & shapes_tmp = shapes.alloc(nb_element*nb_total_points,
					   size_of_shapes,
					   itp_type,
					   ghost_type);

  Array<Real>::matrix_iterator shapes_it =
    shapes_tmp.begin_reinterpret(ElementClass<type>::getNbNodesPerInterpolationElement(), nb_total_points, nb_element);


  Matrix<Real> physical_points_1(spatial_dimension, nb_points_sub_1);
  Matrix<Real> physical_points_2(spatial_dimension, nb_points_sub_2);
  Matrix<Real> parent_natural_coords_1(spatial_dimension, nb_points_sub_1);
  Matrix<Real> parent_natural_coords_2(spatial_dimension, nb_points_sub_2);

  /// intialize the matrices for the parent and subelement coordinates
  UInt nb_nodes_parent_el = ElementClass<parent_type>::getNbNodesPerInterpolationElement();
  UInt nb_nodes_sub_el_1 = ElementClass<sub_type_1>::getNbNodesPerInterpolationElement();
  UInt nb_nodes_sub_el_2 = ElementClass<sub_type_2>::getNbNodesPerInterpolationElement();
  Matrix<Real> parent_coords(spatial_dimension, nb_nodes_parent_el);
  Matrix<Real> sub_el_1_coords(spatial_dimension, nb_nodes_sub_el_1);
  Matrix<Real> sub_el_2_coords(spatial_dimension, nb_nodes_sub_el_2);

  /// loop over all elements of the given type and compute the shape functions
  Vector<Real> all_shapes(size_of_shapes);

  for (UInt elem = 0; elem < nb_element; ++elem, ++shapes_it, ++x_it, ++natural_coords_parent_it) {
    Matrix<Real> & N = *shapes_it;
    const Matrix<Real> & X = *x_it;
    Matrix<Real> & nc_parent = *natural_coords_parent_it;

    /// map the sub element control points into the parent reference domain
    ElementClass<type>::mapFromSubRefToParentRef(X, sub_el_1_coords, parent_coords, sub_1_shapes, physical_points_1, parent_natural_coords_1, 0);
    ElementClass<type>::mapFromSubRefToParentRef(X, sub_el_2_coords, parent_coords, sub_2_shapes, physical_points_2, parent_natural_coords_2, 1);

    /// compute the parent shape functions on all control points
    ElementClass<sub_type_1>::computeShapes(parent_natural_coords_1, parent_1_shapes);
    ElementClass<sub_type_1>::computeShapes(parent_natural_coords_2, parent_2_shapes);

    /// copy the results into the shape functions iterator and natural coords iterator
    for (UInt i = 0; i < nb_points_sub_1; ++i) {
      ElementClass<type>::assembleShapes(parent_1_shapes(i), sub_1_shapes(i), all_shapes, 0);
      N(i) = all_shapes;
      nc_parent(i) = parent_natural_coords_1(i);
    }
    for (UInt i = 0; i < nb_points_sub_2; ++i) {
      ElementClass<type>::assembleShapes(parent_2_shapes(i), sub_2_shapes(i), all_shapes, 1);
      N(i + nb_points_sub_1) = all_shapes;
      
      ///N(i + nb_points_sub_2) = all_shapes;
      nc_parent(i + nb_points_sub_1) = parent_natural_coords_2(i);
    }
      
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_igfem>::precomputeShapeDerivativesOnControlPoints(const Array<Real> & nodes, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// typedef for the two subelement_types and the parent element type
  const ElementType sub_type_1 = ElementClassProperty<type>::sub_element_type_1;
  const ElementType sub_type_2 = ElementClassProperty<type>::sub_element_type_2;
  const ElementType parent_type = ElementClassProperty<type>::parent_element_type;

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  UInt spatial_dimension = mesh.getSpatialDimension();

  /// get the integration points for the subelements 
  Matrix<Real> & natural_coords_sub_1 = control_points(sub_type_1, ghost_type);
  Matrix<Real> & natural_coords_sub_2 = control_points(sub_type_2, ghost_type);
  /// store the number of quadrature points on each subelement and the toal number 
  UInt nb_points_sub_1 = natural_coords_sub_1.cols();
  UInt nb_points_sub_2 = natural_coords_sub_2.cols();
  UInt nb_points_total = nb_points_sub_1 + nb_points_sub_2;

  UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerInterpolationElement();
  UInt size_of_shapesd      = ElementClass<type>::getShapeDerivativesSize();
   
  /// intialize the matrices for the parent and subelement coordinates
  UInt nb_nodes_parent_el = ElementClass<parent_type>::getNbNodesPerInterpolationElement();
  UInt nb_nodes_sub_el_1 = ElementClass<sub_type_1>::getNbNodesPerInterpolationElement();
  UInt nb_nodes_sub_el_2 = ElementClass<sub_type_2>::getNbNodesPerInterpolationElement();
  Matrix<Real> parent_coords(spatial_dimension, nb_nodes_parent_el);
  Matrix<Real> sub_el_1_coords(spatial_dimension, nb_nodes_sub_el_1);
  Matrix<Real> sub_el_2_coords(spatial_dimension, nb_nodes_sub_el_2);

  UInt nb_element = mesh.getConnectivity(type, ghost_type).getSize();
  Array<Real> & shapes_derivatives_tmp = shapes_derivatives.alloc(nb_element*nb_points_total,
  								   size_of_shapesd,
  								   itp_type,
  								   ghost_type);

  /// get an iterator to the coordiantes of the elements
  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el,
				       type, ghost_type);

  Real * shapesd_val = shapes_derivatives_tmp.storage();
  Array<Real>::matrix_iterator x_it = x_el.begin(spatial_dimension,
  								  nb_nodes_per_element);

  /// get an iterator to the control points of the parent element
  Array<Real> & natural_coords_parent = igfem_control_points(type, ghost_type);
  Array<Real>::matrix_iterator natural_coords_parent_it = natural_coords_parent.begin_reinterpret(spatial_dimension, nb_points_total, nb_element);
  
  Tensor3<Real> B_sub_1(spatial_dimension,  nb_nodes_sub_el_1, nb_points_sub_1);
  Tensor3<Real> B_sub_2(spatial_dimension, nb_nodes_sub_el_2, nb_points_sub_2);
  Tensor3<Real> B_parent(spatial_dimension, nb_nodes_parent_el, nb_points_total);
  /// assemble the shape derivatives 
  Matrix<Real> all_shapes(spatial_dimension, nb_nodes_per_element);
  for (UInt elem = 0; elem < nb_element; ++elem, ++x_it, ++natural_coords_parent_it) {
    Matrix<Real> & X = *x_it;
    Matrix<Real> & nc_parent = *natural_coords_parent_it;

    Tensor3<Real> B(shapesd_val,
		    spatial_dimension, nb_nodes_per_element, nb_points_total);
  
    /// get the coordinates of the two sub elements and the parent element
    ElementClass<type>::getSubElementCoords(X, sub_el_1_coords, 0);
    ElementClass<type>::getSubElementCoords(X, sub_el_2_coords, 1);
    ElementClass<type>::getParentCoords(X, parent_coords);

    /// compute the subelements' shape derivatives and the parent shape derivatives
    computeShapeDerivativesOnCPointsByElement<sub_type_1>(sub_el_1_coords,
							  natural_coords_sub_1,
							  B_sub_1);
    computeShapeDerivativesOnCPointsByElement<sub_type_2>(sub_el_2_coords,
							  natural_coords_sub_1,
							  B_sub_2);
    computeShapeDerivativesOnCPointsByElement<parent_type>(parent_coords,
							  nc_parent,
							  B_parent);

    for (UInt i = 0; i < nb_points_sub_1; ++i) {
      ElementClass<type>::assembleShapeDerivatives(B_parent(i), B_sub_1(i), all_shapes, 0);
      B(i) = all_shapes;
    }

    for (UInt i = 0; i < nb_points_sub_2; ++i) {
      ElementClass<type>::assembleShapeDerivatives(B_parent(i), B_sub_2(i), all_shapes, 1);
      B(i + nb_points_sub_1) = all_shapes;
    }
    
    shapesd_val += size_of_shapesd*nb_points_total;
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

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrange<_ek_igfem>::fieldTimesShapes(const Array<Real> & field,
					   Array<Real> & field_times_shapes,
					   GhostType ghost_type) const {
  AKANTU_DEBUG_IN();

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

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */




