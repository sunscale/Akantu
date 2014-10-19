/**
 * @file   fe_engine_template_tmpl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Nov 05 2012
 * @date last modification: Mon Jul 07 2014
 *
 * @brief  Template implementation of FEEngineTemplate
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

#include "aka_common.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
FEEngineTemplate<I, S, kind>::FEEngineTemplate(Mesh & mesh, UInt spatial_dimension,
					       ID id, MemoryID memory_id) :
  FEEngine(mesh,spatial_dimension,id,memory_id),
  integrator(mesh, id, memory_id),
  shape_functions(mesh, id, memory_id) { }

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
FEEngineTemplate<I, S, kind>::~FEEngineTemplate() { }

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct GradientOnQuadraturePointsHelper {
  template <class S>
  static void call(const S & shape_functions,	
		   Mesh & mesh,
		   const Array<Real> & u,
		   Array<Real> & nablauq,
		   const UInt nb_degree_of_freedom,
		   const ElementType & type,
		   const GhostType & ghost_type,
		   const Array<UInt> & filter_elements) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

#define COMPUTE_GRADIENT(type)						\
  if (element_dimension == ElementClass<type>::getSpatialDimension())	\
    shape_functions.template gradientOnControlPoints<type>(u,		\
							   nablauq,	\
							   nb_degree_of_freedom, \
							   ghost_type,	\
							   filter_elements);

#define AKANTU_SPECIALIZE_GRADIENT_ON_QUADRATURE_POINTS_HELPER(kind)	\
  template<>								\
  struct GradientOnQuadraturePointsHelper<kind> {			\
    template <class S>							\
    static void call(const S & shape_functions,				\
		     Mesh & mesh,					\
		     const Array<Real> & u,				\
		     Array<Real> & nablauq,				\
		     const UInt nb_degree_of_freedom,			\
		     const ElementType & type,				\
		     const GhostType & ghost_type,			\
		     const Array<UInt> & filter_elements) {		\
      UInt element_dimension = mesh.getSpatialDimension(type);		\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(COMPUTE_GRADIENT, kind);		\
    }									\
  };
#define INTEREST_LIST AKANTU_GENERATE_KIND_LIST(AKANTU_REGULAR_KIND AKANTU_COHESIVE_KIND AKANTU_IGFEM_KIND)

AKANTU_BOOST_ALL_KIND_LIST(AKANTU_SPECIALIZE_GRADIENT_ON_QUADRATURE_POINTS_HELPER, \
			   INTEREST_LIST)

#undef AKANTU_SPECIALIZE_GRADIENT_ON_QUADRATURE_POINTS_HELPER
#undef COMPUTE_GRADIENT
#undef INTEREST_LIST

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEEngineTemplate<I, S, kind>::gradientOnQuadraturePoints(const Array<Real> &u,
							      Array<Real> &nablauq,
							      const UInt nb_degree_of_freedom,
							      const ElementType & type,
							      const GhostType & ghost_type,
							      const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  if(filter_elements != empty_filter) nb_element = filter_elements.getSize();
  UInt nb_points         = shape_functions.getControlPoints(type, ghost_type).cols();

#ifndef AKANTU_NDEBUG

  UInt element_dimension = mesh.getSpatialDimension(type);

  AKANTU_DEBUG_ASSERT(u.getSize() == mesh.getNbNodes(),
		      "The vector u(" << u.getID()
		      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(u.getNbComponent() == nb_degree_of_freedom ,
		      "The vector u(" << u.getID()
		      << ") has not the good number of component.");

  AKANTU_DEBUG_ASSERT(nablauq.getNbComponent()
		      == nb_degree_of_freedom * element_dimension,
		      "The vector nablauq(" << nablauq.getID()
		      << ") has not the good number of component.");

  // AKANTU_DEBUG_ASSERT(nablauq.getSize() == nb_element * nb_points,
  //	   	      "The vector nablauq(" << nablauq.getID()
  //	   	      << ") has not the good size.");
#endif

  nablauq.resize(nb_element * nb_points);

  GradientOnQuadraturePointsHelper<kind>::call(shape_functions,
					       mesh,
					       u,
					       nablauq,
					       nb_degree_of_freedom,
					       type,
					       ghost_type,
					       filter_elements);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEEngineTemplate<I, S, kind>::initShapeFunctions(const GhostType & ghost_type) {
  initShapeFunctions(mesh.getNodes(), ghost_type);
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEEngineTemplate<I, S, kind>::initShapeFunctions(const Array<Real> & nodes,
						      const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it  = mesh.firstType(element_dimension, ghost_type, kind);
  Mesh::type_iterator end = mesh.lastType(element_dimension, ghost_type, kind);
  for(; it != end; ++it) {
    ElementType type = *it;
    integrator.initIntegrator(nodes, type, ghost_type);
    const Matrix<Real> & control_points =
      getQuadraturePoints(type, ghost_type);
    shape_functions.initShapeFunctions(nodes, control_points, type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct IntegrateHelper { };

#define INTEGRATE(type)						\
  integrator.template integrate<type>(f,			\
				      intf,			\
				      nb_degree_of_freedom,	\
				      ghost_type,		\
				      filter_elements);		\

#define AKANTU_SPECIALIZE_INTEGRATE_HELPER(kind)		\
  template<>							\
  struct IntegrateHelper<kind> {				\
    template <class I>						\
    static void call(const I & integrator,			\
		     const Array<Real> & f,			\
		     Array<Real> &intf,				\
		     UInt nb_degree_of_freedom,			\
		     const ElementType & type,			\
		     const GhostType & ghost_type,		\
		     const Array<UInt> & filter_elements) {	\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(INTEGRATE, kind);	\
    }								\
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_INTEGRATE_HELPER)

#undef AKANTU_SPECIALIZE_INTEGRATE_HELPER
#undef INTEGRATE

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEEngineTemplate<I, S, kind>::integrate(const Array<Real> & f,
					     Array<Real> &intf,
					     UInt nb_degree_of_freedom,
					     const ElementType & type,
					     const GhostType & ghost_type,
					     const Array<UInt> & filter_elements) const{

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  if(filter_elements != empty_filter) nb_element = filter_elements.getSize();
#ifndef AKANTU_NDEBUG

  UInt nb_quadrature_points  = getNbQuadraturePoints(type);

  AKANTU_DEBUG_ASSERT(f.getSize() == nb_element * nb_quadrature_points,
		      "The vector f(" << f.getID() << " size " << f.getSize()
		      << ") has not the good size (" << nb_element << ").");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degree_of_freedom ,
		      "The vector f(" << f.getID()
		      << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(intf.getNbComponent() == nb_degree_of_freedom,
		      "The vector intf(" << intf.getID()
		      << ") has not the good number of component.");
#endif

  intf.resize(nb_element);

  IntegrateHelper<kind>::call(integrator,
			      f,
			      intf,
			      nb_degree_of_freedom,
			      type,
			      ghost_type,
			      filter_elements);

}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct IntegrateScalarHelper { };

#define INTEGRATE(type)							\
  integral = integrator.template integrate<type>(f,			\
						 ghost_type, filter_elements);

#define AKANTU_SPECIALIZE_INTEGRATE_SCALAR_HELPER(kind)		\
  template<>							\
  struct IntegrateScalarHelper<kind> {				\
    template <class I>						\
    static Real call(const I & integrator,			\
		     const Array<Real> & f,			\
		     const ElementType & type,			\
		     const GhostType & ghost_type,		\
		     const Array<UInt> & filter_elements) {	\
      Real integral = 0.;					\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(INTEGRATE, kind);	\
      return integral;						\
    }								\
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_INTEGRATE_SCALAR_HELPER)

#undef AKANTU_SPECIALIZE_INTEGRATE_SCALAR_HELPER
#undef INTEGRATE

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
Real FEEngineTemplate<I, S, kind>::integrate(const Array<Real> & f,
					     const ElementType & type,
					     const GhostType & ghost_type,
					     const Array<UInt> & filter_elements) const{
  AKANTU_DEBUG_IN();

#ifndef AKANTU_NDEBUG
  //   std::stringstream sstr; sstr << ghost_type;
  //   AKANTU_DEBUG_ASSERT(sstr.str() == nablauq.getTag(),
  //		      "The vector " << nablauq.getID() << " is not taged " << ghost_type);
  UInt nb_element = mesh.getNbElement(type, ghost_type);
  if(filter_elements != empty_filter) nb_element = filter_elements.getSize();

  UInt nb_quadrature_points  = getNbQuadraturePoints(type, ghost_type);

  AKANTU_DEBUG_ASSERT(f.getSize() == nb_element * nb_quadrature_points,
		      "The vector f(" << f.getID()
		      << ") has not the good size. (" << f.getSize() << "!=" << nb_quadrature_points * nb_element << ")");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == 1,
		      "The vector f(" << f.getID()
		      << ") has not the good number of component.");
#endif

  Real integral = IntegrateScalarHelper<kind>::call(integrator, f, type, ghost_type, filter_elements);
  AKANTU_DEBUG_OUT();
  return integral;
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct IntegrateScalarOnOneElementHelper { };

#define INTEGRATE(type)						\
  res = integrator.template integrate<type>(f,			\
					    index, ghost_type);

#define AKANTU_SPECIALIZE_INTEGRATE_SCALAR_ON_ONE_ELEMENT_HELPER(kind)	\
  template<>								\
  struct IntegrateScalarOnOneElementHelper<kind> {			\
    template <class I>							\
    static Real call(const I & integrator,				\
		     const Vector<Real> & f,				\
		     const ElementType & type,				\
		     UInt index,					\
		     const GhostType & ghost_type) {			\
      Real res = 0.;							\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(INTEGRATE, kind);		\
      return res;							\
    }									\
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_INTEGRATE_SCALAR_ON_ONE_ELEMENT_HELPER)

#undef AKANTU_SPECIALIZE_INTEGRATE_SCALAR_ON_ONE_ELEMENT_HELPER
#undef INTEGRATE

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
Real FEEngineTemplate<I, S, kind>::integrate(const Vector<Real> & f,
					     const ElementType & type,
					     UInt index,
					     const GhostType & ghost_type) const{


  Real res = IntegrateScalarOnOneElementHelper<kind>::call(integrator, f, type, index, ghost_type);
  return res;
}


/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct IntegrateOnQuadraturePointsHelper { };

#define INTEGRATE(type)							\
  integrator.template integrateOnQuadraturePoints<type>(f, intf, nb_degree_of_freedom, \
							ghost_type, filter_elements);

#define AKANTU_SPECIALIZE_INTEGRATE_ON_QUADRATURE_POINTS_HELPER(kind)	\
  template<>								\
  struct IntegrateOnQuadraturePointsHelper<kind> {			\
    template <class I>							\
    static void call(const I & integrator,				\
		     const Array<Real> & f,				\
		     Array<Real> & intf,				\
		     UInt nb_degree_of_freedom,				\
		     const ElementType & type,				\
		     const GhostType & ghost_type,			\
		     const Array<UInt> & filter_elements) {		\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(INTEGRATE, kind);		\
    }									\
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_INTEGRATE_ON_QUADRATURE_POINTS_HELPER)

#undef AKANTU_SPECIALIZE_INTEGRATE_ON_QUADRATURE_POINTS_HELPER
#undef INTEGRATE

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEEngineTemplate<I, S, kind>::integrateOnQuadraturePoints(const Array<Real> & f,
							       Array<Real> &intf,
							       UInt nb_degree_of_freedom,
							       const ElementType & type,
							       const GhostType & ghost_type,
							       const Array<UInt> & filter_elements) const{

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  if(filter_elements != empty_filter) nb_element = filter_elements.getSize();
  UInt nb_quadrature_points  = getNbQuadraturePoints(type);
#ifndef AKANTU_NDEBUG
  //   std::stringstream sstr; sstr << ghost_type;
  //   AKANTU_DEBUG_ASSERT(sstr.str() == nablauq.getTag(),
  //		      "The vector " << nablauq.getID() << " is not taged " << ghost_type);


  AKANTU_DEBUG_ASSERT(f.getSize() == nb_element * nb_quadrature_points,
		      "The vector f(" << f.getID() << " size " << f.getSize()
		      << ") has not the good size (" << nb_element << ").");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degree_of_freedom ,
		      "The vector f(" << f.getID()
		      << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(intf.getNbComponent() == nb_degree_of_freedom,
		      "The vector intf(" << intf.getID()
		      << ") has not the good number of component.");
#endif

  intf.resize(nb_element*nb_quadrature_points);
  IntegrateOnQuadraturePointsHelper<kind>::call(integrator, f, intf, nb_degree_of_freedom, type, ghost_type, filter_elements);
}


/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct InterpolateOnQuadraturePointsHelper {
  template <class S>
  static void call(const S & shape_functions,
		   const Array<Real> &u,
		   Array<Real> &uq,
		   const UInt nb_degree_of_freedom,
		   const ElementType & type,
		   const GhostType & ghost_type,
		   const Array<UInt> & filter_elements) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};


#define INTERPOLATE(type)						\
  shape_functions.template interpolateOnControlPoints<type>(u,		\
							    uq,		\
							    nb_degree_of_freedom, \
							    ghost_type, \
							    filter_elements);

#define AKANTU_SPECIALIZE_INTERPOLATE_ON_QUADRATURE_POINTS_HELPER(kind)	\
  template<>								\
  struct InterpolateOnQuadraturePointsHelper<kind> {			\
    template <class S>							\
    static void call(const S & shape_functions,				\
		     const Array<Real> & u,				\
		     Array<Real> & uq,					\
		     const UInt nb_degree_of_freedom,			\
		     const ElementType & type,				\
		     const GhostType & ghost_type,			\
		     const Array<UInt> & filter_elements) {		\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(INTERPOLATE, kind);		\
    }									\
  };
#define INTEREST_LIST AKANTU_GENERATE_KIND_LIST(AKANTU_REGULAR_KIND AKANTU_COHESIVE_KIND AKANTU_IGFEM_KIND)

AKANTU_BOOST_ALL_KIND_LIST(AKANTU_SPECIALIZE_INTERPOLATE_ON_QUADRATURE_POINTS_HELPER, \
			   INTEREST_LIST)

#undef AKANTU_SPECIALIZE_INTERPOLATE_ON_QUADRATURE_POINTS_HELPER
#undef INTERPOLATE
#undef INTEREST_LIST

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEEngineTemplate<I, S, kind>::interpolateOnQuadraturePoints(const Array<Real> &u,
								 Array<Real> &uq,
								 const UInt nb_degree_of_freedom,
								 const ElementType & type,
								 const GhostType & ghost_type,
								 const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_points = shape_functions.getControlPoints(type, ghost_type).cols();
  UInt nb_element = mesh.getNbElement(type, ghost_type);
  if(filter_elements != empty_filter) nb_element = filter_elements.getSize();
#ifndef AKANTU_NDEBUG

  AKANTU_DEBUG_ASSERT(u.getSize() == mesh.getNbNodes(),
		      "The vector u(" << u.getID()
		      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(u.getNbComponent() == nb_degree_of_freedom ,
		      "The vector u(" << u.getID()
		      << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(uq.getNbComponent() == nb_degree_of_freedom,
		      "The vector uq(" << uq.getID()
		      << ") has not the good number of component.");
#endif

  uq.resize(nb_element * nb_points);


  InterpolateOnQuadraturePointsHelper<kind>::call(shape_functions,
						  u,
						  uq,
						  nb_degree_of_freedom,
						  type,
						  ghost_type,
						  filter_elements);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEEngineTemplate<I, S, kind>::interpolateOnQuadraturePoints(const Array<Real> & u,
								 ElementTypeMapArray<Real> & uq,
								 const ElementTypeMapArray<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  const Array<UInt> * filter = NULL;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType ghost_type = *gt;
    ElementTypeMapArray<Real>::type_iterator it   = uq.firstType(_all_dimensions, ghost_type, kind);
    ElementTypeMapArray<Real>::type_iterator last = uq.lastType(_all_dimensions, ghost_type, kind);

    for (; it != last; ++it) {
      ElementType type = *it;

      UInt nb_quad_per_element = getNbQuadraturePoints(type, ghost_type);

      UInt nb_element = 0;

      if (filter_elements) {
	filter = &((*filter_elements)(type, ghost_type));
	nb_element = filter->getSize();
      }
      else {
	filter = &empty_filter;
	nb_element = mesh.getNbElement(type, ghost_type);
      }

      UInt nb_tot_quad = nb_quad_per_element * nb_element;

      Array<Real> & quad = uq(type, ghost_type);
      quad.resize(nb_tot_quad);

      interpolateOnQuadraturePoints(u,
				    quad,
				    mesh.getSpatialDimension(),
				    type,
				    ghost_type,
				    *filter);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEEngineTemplate<I, S, kind>::computeNormalsOnControlPoints(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  computeNormalsOnControlPoints(mesh.getNodes(),
				ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEEngineTemplate<I, S, kind>::computeNormalsOnControlPoints(const Array<Real> & field,
								 const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  //  Real * coord = mesh.getNodes().storage();
  UInt spatial_dimension = mesh.getSpatialDimension();

  //allocate the normal arrays
  Mesh::type_iterator it  = mesh.firstType(element_dimension, ghost_type, kind);
  Mesh::type_iterator end = mesh.lastType(element_dimension, ghost_type, kind);
  for(; it != end; ++it) {
    ElementType type = *it;
    UInt size = mesh.getNbElement(type, ghost_type);
    if(normals_on_quad_points.exists(type, ghost_type)) {
      normals_on_quad_points(type, ghost_type).resize(size);
    } else {
      normals_on_quad_points.alloc(size, spatial_dimension, type, ghost_type);
    }
  }

  //loop over the type to build the normals
  it  = mesh.firstType(element_dimension, ghost_type, kind);
  for(; it != end; ++it) {
    Array<Real> & normals_on_quad = normals_on_quad_points(*it, ghost_type);
    computeNormalsOnControlPoints(field, normals_on_quad, *it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct ComputeNormalsOnControlPoints {
  template <template <ElementKind> class I,
	    template <ElementKind> class S,
	    ElementKind k>
  static void call(const FEEngineTemplate<I, S, k> & fem,
		   const Array<Real> & field,
		   Array<Real> & normal,
		   const ElementType & type,
		   const GhostType & ghost_type) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

#define COMPUTE_NORMALS_ON_QUAD(type)					\
  fem.template computeNormalsOnControlPoints<type>(field, normal, ghost_type);

#define AKANTU_SPECIALIZE_COMPUTE_NORMALS_ON_CONTROL_POINTS(kind)	\
  template<>								\
  struct ComputeNormalsOnControlPoints<kind> {				\
    template <template <ElementKind> class I,				\
	      template <ElementKind> class S,				\
	      ElementKind k>						\
    static void call(const FEEngineTemplate<I, S, k> & fem,		\
		     const Array<Real> & field,				\
		     Array<Real> & normal,				\
		     const ElementType & type,				\
		     const GhostType & ghost_type) {			\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(COMPUTE_NORMALS_ON_QUAD, kind);	\
    }									\
    };

#define INTEREST_LIST AKANTU_GENERATE_KIND_LIST(AKANTU_REGULAR_KIND AKANTU_COHESIVE_KIND AKANTU_IGFEM_KIND)

AKANTU_BOOST_ALL_KIND_LIST(AKANTU_SPECIALIZE_COMPUTE_NORMALS_ON_CONTROL_POINTS, \
			   INTEREST_LIST)

#undef AKANTU_SPECIALIZE_COMPUTE_NORMALS_ON_CONTROL_POINTS
#undef COMPUTE_NORMALS_ON_QUAD
#undef INTEREST_LIST

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEEngineTemplate<I, S, kind>::computeNormalsOnControlPoints(const Array<Real> & field,
								 Array<Real> & normal,
								 const ElementType & type,
								 const GhostType & ghost_type) const {
  ComputeNormalsOnControlPoints<kind>::call(*this,
					    field,
					    normal,
					    type,
					    ghost_type);
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
template<ElementType type>
void FEEngineTemplate<I, S, kind>::computeNormalsOnControlPoints(const Array<Real> & field,
								 Array<Real> & normal,
								 const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_nodes_per_element  = Mesh::getNbNodesPerElement(type);
  UInt nb_points = getNbQuadraturePoints(type, ghost_type);

  UInt nb_element = mesh.getConnectivity(type, ghost_type).getSize();
  normal.resize(nb_element * nb_points);
  Array<Real>::matrix_iterator normals_on_quad = normal.begin_reinterpret(spatial_dimension,
									  nb_points,
									  nb_element);
  Array<Real> f_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, field, f_el, type, ghost_type);

  const Matrix<Real> & quads =
    integrator. template getQuadraturePoints<type>(ghost_type);

  Array<Real>::matrix_iterator f_it = f_el.begin(spatial_dimension, nb_nodes_per_element);

  for (UInt elem = 0; elem < nb_element; ++elem) {
    ElementClass<type>::computeNormalsOnNaturalCoordinates(quads,
							   *f_it,
							   *normals_on_quad);
    ++normals_on_quad;
    ++f_it;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/* Matrix lumping functions                                                   */
/* -------------------------------------------------------------------------- */

/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct AssembleLumpedTemplateHelper { };

#define ASSEMBLE_LUMPED(type)                                           \
  fem.template assembleLumpedTemplate<type>(field_1, nb_degree_of_freedom,lumped, equation_number,ghost_type)

#define AKANTU_SPECIALIZE_ASSEMBLE_HELPER(kind)			\
  template<>							\
  struct AssembleLumpedTemplateHelper<kind> {			\
    template <template <ElementKind> class I,			\
	      template <ElementKind> class S,			\
	      ElementKind k>					\
    static void call(const FEEngineTemplate<I, S, k> & fem,	\
		     const Array<Real> & field_1,		\
		     UInt nb_degree_of_freedom,			\
		     Array<Real> & lumped,			\
		     const Array<Int> & equation_number,	\
		     ElementType type,				\
		     const GhostType & ghost_type) {		\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(ASSEMBLE_LUMPED, kind);	\
    }								\
    };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_ASSEMBLE_HELPER)

#undef AKANTU_SPECIALIZE_ASSEMBLE_HELPER
#undef ASSEMBLE_LUMPED

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEEngineTemplate<I, S, kind>::assembleFieldLumped(const Array<Real> & field_1,
						       UInt nb_degree_of_freedom,
						       Array<Real> & lumped,
						       const Array<Int> & equation_number,
						       ElementType type,
						       const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  AssembleLumpedTemplateHelper<kind>::call(*this, field_1,
					   nb_degree_of_freedom,
					   lumped,
					   equation_number,
					   type,
					   ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct AssembleFieldMatrixHelper { };

#define ASSEMBLE_MATRIX(type)                                           \
  fem.template assembleFieldMatrix<type>(field_1, nb_degree_of_freedom,	\
					 matrix,			\
					 ghost_type)

#define AKANTU_SPECIALIZE_ASSEMBLE_FIELD_MATRIX_HELPER(kind)	\
  template<>							\
  struct AssembleFieldMatrixHelper<kind> {			\
    template <template <ElementKind> class I,			\
	      template <ElementKind> class S,			\
	      ElementKind k>					\
    static void call(const FEEngineTemplate<I, S, k> & fem,	\
		     const Array<Real> & field_1,		\
		     UInt nb_degree_of_freedom,			\
		     SparseMatrix & matrix,			\
		     ElementType type,				\
		     const GhostType & ghost_type) {		\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(ASSEMBLE_MATRIX, kind);	\
    }								\
    };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_ASSEMBLE_FIELD_MATRIX_HELPER)

#undef AKANTU_SPECIALIZE_ASSEMBLE_FIELD_MATRIX_HELPER
#undef ASSEMBLE_MATRIX

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
void FEEngineTemplate<I, S, kind>::assembleFieldMatrix(const Array<Real> & field_1,
						       UInt nb_degree_of_freedom,
						       SparseMatrix & matrix,
						       ElementType type,
						       const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();
  AssembleFieldMatrixHelper<kind>::call(*this, field_1, nb_degree_of_freedom, matrix, type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
template <ElementType type>
void FEEngineTemplate<I, S, kind>::assembleLumpedTemplate(const Array<Real> & field_1,
							  UInt nb_degree_of_freedom,
							  Array<Real> & lumped,
							  const Array<Int> & equation_number,
							  const GhostType & ghost_type) const {
  this->template assembleLumpedRowSum<type>(field_1, nb_degree_of_freedom,lumped, equation_number,ghost_type);
}

/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j dV = \int \rho \varphi_i dV @f$
 */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
template <ElementType type>
void FEEngineTemplate<I, S, kind>::assembleLumpedRowSum(const Array<Real> & field_1,
							UInt nb_degree_of_freedom,
							Array<Real> & lumped,
							const Array<Int> & equation_number,
							const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt shapes_size = ElementClass<type>::getShapeSize();

  Array<Real> * field_times_shapes = new Array<Real>(0, shapes_size * nb_degree_of_freedom);
  Array<Real> * field = new Array<Real>(field_1.getSize(), nb_degree_of_freedom);

  Array<Real>::const_scalar_iterator f1_it  = field_1.begin();
  Array<Real>::const_scalar_iterator f1_end = field_1.end();
  Array<Real>::vector_iterator f_it = field->begin(nb_degree_of_freedom);

  for(;f1_it != f1_end; ++f1_it, ++f_it) {
    f_it->set(*f1_it);
  }

  shape_functions.template fieldTimesShapes<type>(*field, *field_times_shapes, ghost_type);
  delete field;

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  Array<Real> * int_field_times_shapes = new Array<Real>(nb_element, shapes_size * nb_degree_of_freedom,
							 "inte_rho_x_shapes");
  integrator.template integrate<type>(*field_times_shapes, *int_field_times_shapes,
				      nb_degree_of_freedom * shapes_size, ghost_type, empty_filter);
  delete field_times_shapes;

  assembleArray(*int_field_times_shapes, lumped, equation_number,nb_degree_of_freedom, type, ghost_type);
  delete int_field_times_shapes;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = c * M_{ii} = \int_{V_e} \rho dV @f$
 */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
template <ElementType type>
void FEEngineTemplate<I, S, kind>::assembleLumpedDiagonalScaling(const Array<Real> & field_1,
								 UInt nb_degree_of_freedom,
								 Array<Real> & lumped,
								 const Array<Int> & equation_number,
								 const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  const ElementType & type_p1  = ElementClass<type>::getP1ElementType();
  UInt nb_nodes_per_element_p1 = Mesh::getNbNodesPerElement(type_p1);
  UInt nb_nodes_per_element    = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points    = integrator.template getQuadraturePoints<type>(ghost_type).cols();

  UInt nb_element = field_1.getSize() / nb_quadrature_points;

  Real corner_factor = 0;
  Real mid_factor    = 0;

  if(type == _triangle_6) {
    corner_factor = 1./12.;
    mid_factor    = 1./4.;
  }

  if (type == _tetrahedron_10) {
    corner_factor = 1./32.;
    mid_factor    = 7./48.;
  }

  if (type == _quadrangle_8) {
    corner_factor = 1./36.;
    mid_factor    = 8./36.;
  }

  if (nb_element == 0) {
    AKANTU_DEBUG_OUT();
    return;
  }

  /// compute @f$ \int \rho dV = \rho V @f$ for each element
  Array<Real> * int_field_1 = new Array<Real>(field_1.getSize(), 1,
					      "inte_rho_x_1");
  integrator.template integrate<type>(field_1, *int_field_1, 1, ghost_type, empty_filter);

  /// distribute the mass of the element to the nodes
  Array<Real> * lumped_per_node = new Array<Real>(nb_element, nb_degree_of_freedom * nb_nodes_per_element, "mass_per_node");
  Array<Real>::const_scalar_iterator int_field_1_it = int_field_1->begin();
  Array<Real>::matrix_iterator lumped_per_node_it
    = lumped_per_node->begin(nb_degree_of_freedom, nb_nodes_per_element);

  for (UInt e = 0; e < nb_element; ++e) {
    Real lmass = *int_field_1_it * corner_factor;
    for (UInt n = 0; n < nb_nodes_per_element_p1; ++n) {
      Vector<Real> l = (*lumped_per_node_it)(n);
      l.set(lmass); /// corner points
    }

    lmass = *int_field_1_it * mid_factor;
    for (UInt n = nb_nodes_per_element_p1; n < nb_nodes_per_element; ++n) {
      Vector<Real> l = (*lumped_per_node_it)(n);
      l.set(lmass); /// mid points
    }

    ++int_field_1_it;
    ++lumped_per_node_it;
  }
  delete int_field_1;

  //  lumped_per_node->extendComponentsInterlaced(nb_degree_of_freedom,1);
  assembleArray(*lumped_per_node, lumped, equation_number, nb_degree_of_freedom, type, ghost_type);
  delete lumped_per_node;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j dV = \int \rho \varphi_i dV @f$
 */
template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
template <ElementType type>
void FEEngineTemplate<I, S, kind>::assembleFieldMatrix(const Array<Real> & field_1,
						       UInt nb_degree_of_freedom,
						       SparseMatrix & matrix,
						       const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt vect_size   = field_1.getSize();
  UInt shapes_size = ElementClass<type>::getShapeSize();
  UInt lmat_size   = nb_degree_of_freedom * shapes_size;

  const Array<Real> & shapes = shape_functions.getShapes(type,ghost_type);
  Array<Real> * modified_shapes = new Array<Real>(vect_size, lmat_size * nb_degree_of_freedom);
  modified_shapes->clear();
  Array<Real> * local_mat = new Array<Real>(vect_size, lmat_size * lmat_size);

  Array<Real>::matrix_iterator shape_vect  = modified_shapes->begin(nb_degree_of_freedom, lmat_size);
  Real * sh  = shapes.storage();
  for(UInt q = 0; q < vect_size; ++q) {
    Real * msh = shape_vect->storage();
    for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
      Real * msh_tmp = msh + d * (lmat_size + 1);
      for (UInt s = 0; s < shapes_size; ++s) {
	*msh_tmp = sh[s];
	msh_tmp += nb_degree_of_freedom;
      }
    }
    ++shape_vect;
    sh += shapes_size;
  }

  shape_vect  = modified_shapes->begin(nb_degree_of_freedom, lmat_size);
  Array<Real>::matrix_iterator lmat = local_mat->begin(lmat_size, lmat_size);
  Real * field_val = field_1.storage();

  for(UInt q = 0; q < vect_size; ++q) {
    (*lmat).mul<true, false>(*shape_vect, *shape_vect, *field_val);
    ++lmat; ++shape_vect; ++field_val;
  }

  delete modified_shapes;

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  Array<Real> * int_field_times_shapes = new Array<Real>(nb_element, lmat_size * lmat_size,
							 "inte_rho_x_shapes");
  integrator.template integrate<type>(*local_mat, *int_field_times_shapes,
				      lmat_size * lmat_size, ghost_type, empty_filter);
  delete local_mat;

  assembleMatrix(*int_field_times_shapes, matrix, nb_degree_of_freedom, type, ghost_type);
  delete int_field_times_shapes;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct InverseMapHelper {
  template <class S>
  static void call(const S & shape_functions,
		   const Vector<Real> & real_coords,
		   UInt element,
		   const ElementType & type,
		   Vector<Real> & natural_coords,
		   const GhostType & ghost_type) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

#define INVERSE_MAP(type)						\
  shape_functions.template inverseMap<type>(real_coords, element, natural_coords, ghost_type); \

#define AKANTU_SPECIALIZE_INVERSE_MAP_HELPER(kind)		\
  template<>							\
  struct InverseMapHelper<kind> {				\
    template <class S>						\
    static void call(const S & shape_functions,			\
		     const Vector<Real> & real_coords,		\
		     UInt element,				\
		     const ElementType & type,			\
		     Vector<Real> & natural_coords,		\
		     const GhostType & ghost_type) {		\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(INVERSE_MAP, kind);	\
    }								\
  };
#define INTEREST_LIST AKANTU_GENERATE_KIND_LIST(AKANTU_REGULAR_KIND AKANTU_COHESIVE_KIND AKANTU_IGFEM_KIND)

AKANTU_BOOST_ALL_KIND_LIST(AKANTU_SPECIALIZE_INVERSE_MAP_HELPER, \
			   INTEREST_LIST)

#undef AKANTU_SPECIALIZE_INVERSE_MAP_HELPER
#undef INVERSE_MAP
#undef INTEREST_LIST

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
inline void FEEngineTemplate<I, S, kind>::inverseMap(const Vector<Real> & real_coords,
						     UInt element,
						     const ElementType & type,
						     Vector<Real> & natural_coords,
						     const GhostType & ghost_type) const{

  AKANTU_DEBUG_IN();

  InverseMapHelper<kind>::call(shape_functions, real_coords, element, type, natural_coords, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct ContainsHelper {
  template <class S>
  static void call(const S & shape_functions,	
		   const Vector<Real> & real_coords,
		   UInt element,
		   const ElementType & type,
		   const GhostType & ghost_type) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

#define CONTAINS(type)							\
  contain = shape_functions.template contains<type>(real_coords, element, ghost_type); \

#define AKANTU_SPECIALIZE_CONTAINS_HELPER(kind)		\
  template<>						\
  struct ContainsHelper<kind> {				\
    template <template <ElementKind> class S,		\
	      ElementKind k>				\
    static bool call(const S<k> & shape_functions,	\
		     const Vector<Real> & real_coords,	\
		     UInt element,			\
		     const ElementType & type,		\
		     const GhostType & ghost_type) {	\
      bool contain = false;				\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(CONTAINS, kind);	\
      return contain;					\
    }							\
    };
#define INTEREST_LIST AKANTU_GENERATE_KIND_LIST(AKANTU_REGULAR_KIND AKANTU_COHESIVE_KIND AKANTU_IGFEM_KIND)

AKANTU_BOOST_ALL_KIND_LIST(AKANTU_SPECIALIZE_CONTAINS_HELPER, \
			   INTEREST_LIST)

#undef AKANTU_SPECIALIZE_CONTAINS_HELPER
#undef CONTAINS
#undef INTEREST_LIST

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
inline bool FEEngineTemplate<I, S, kind>::contains(const Vector<Real> & real_coords,
						   UInt element,
						   const ElementType & type,
						   const GhostType & ghost_type) const{
  return ContainsHelper<kind>::call(shape_functions, real_coords, element, type, ghost_type);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct ComputeShapesHelper { };

#define COMPUTE_SHAPES(type)						\
  shape_functions.template computeShapes<type>(real_coords,element,shapes,ghost_type); \

#define AKANTU_SPECIALIZE_COMPUTE_SHAPES_HELPER(kind)		\
  template<>							\
  struct ComputeShapesHelper<kind> {				\
    template <class S>						\
    static void call(const S & shape_functions,			\
		     const Vector<Real> & real_coords,		\
		     UInt element,				\
		     const ElementType type,			\
		     Vector<Real> & shapes,			\
		     const GhostType & ghost_type) {		\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(COMPUTE_SHAPES, kind);	\
    }								\
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_COMPUTE_SHAPES_HELPER)

#undef AKANTU_SPECIALIZE_ASSEMBLE_COMPUTE_SHAPES_HELPER
#undef COMPUTE_SHAPES

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
inline void FEEngineTemplate<I, S, kind>::computeShapes(const Vector<Real> & real_coords,
							UInt element,
							const ElementType & type,
							Vector<Real> & shapes,
							const GhostType & ghost_type) const{

  AKANTU_DEBUG_IN();
  
  ComputeShapesHelper<kind>::call(shape_functions, real_coords, element, type, shapes, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct GetNbQuadraturePointsHelper { };

#define GET_NB_QUAD(type)						\
  nb_quad_points =							\
    integrator. template getQuadraturePoints<type>(ghost_type).cols();

#define AKANTU_SPECIALIZE_GET_NB_QUADRATURE_POINTS_HELPER(kind)	\
  template<>							\
  struct GetNbQuadraturePointsHelper<kind> {			\
    template <template <ElementKind> class I,			\
	      ElementKind k>					\
    static UInt call(const I<k> & integrator,			\
		     const ElementType type,			\
		     const GhostType & ghost_type) {		\
      UInt nb_quad_points = 0;					\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(GET_NB_QUAD, kind);	\
      return nb_quad_points;					\
    }								\
    };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_GET_NB_QUADRATURE_POINTS_HELPER)

#undef AKANTU_SPECIALIZE_GET_NB_QUADRATURE_POINTS_HELPER
#undef GET_NB_QUAD

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
inline UInt FEEngineTemplate<I, S, kind>::getNbQuadraturePoints(const ElementType & type,
								const GhostType & ghost_type) const {
  return GetNbQuadraturePointsHelper<kind>::call(integrator, type, ghost_type);
}


/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct GetShapesHelper { };

#define GET_SHAPES(type)                                \
  ret = &(shape_functions.getShapes(type, ghost_type));

#define AKANTU_SPECIALIZE_GET_SHAPES_HELPER(kind)		\
  template<>							\
  struct GetShapesHelper<kind> {				\
  template <class S>						\
  static const Array<Real> & call(const S& shape_functions,	\
			    const ElementType type,		\
			    const GhostType & ghost_type) {	\
    const Array<Real> * ret = NULL;				\
    AKANTU_BOOST_KIND_ELEMENT_SWITCH(GET_SHAPES, kind);		\
    return *ret;						\
  }								\
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_GET_SHAPES_HELPER)

#undef AKANTU_SPECIALIZE_GET_SHAPES_HELPER
#undef GET_SHAPES

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
inline const Array<Real> & FEEngineTemplate<I, S, kind>::getShapes(const ElementType & type,
								   const GhostType & ghost_type,
								   __attribute__((unused)) UInt id) const {
  return GetShapesHelper<kind>::call(shape_functions, type, ghost_type);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct GetShapesDerivativesHelper {
  template <template <ElementKind> class S,
	    ElementKind k>
  static const Array<Real> & call(const S<k> & shape_functions,
				  const ElementType & type,
				  const GhostType & ghost_type,
				  UInt id) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

#define GET_SHAPES_DERIVATIVES(type)					\
  ret = &(shape_functions.getShapesDerivatives(type, ghost_type));

#define AKANTU_SPECIALIZE_GET_SHAPES_DERIVATIVES_HELPER(kind)		\
  template<>								\
  struct GetShapesDerivativesHelper<kind> {				\
    template <template <ElementKind> class S,				\
	      ElementKind k>						\
    static const Array<Real> & call(const S<k> & shape_functions,	\
				    const ElementType type,		\
				    const GhostType & ghost_type,	\
				    UInt id) {				\
      const Array<Real> * ret = NULL;					\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(GET_SHAPES_DERIVATIVES, kind);	\
      return *ret;							\
    }									\
    };

#define INTEREST_LIST AKANTU_GENERATE_KIND_LIST(AKANTU_REGULAR_KIND AKANTU_COHESIVE_KIND AKANTU_IGFEM_KIND)

AKANTU_BOOST_ALL_KIND_LIST(AKANTU_SPECIALIZE_GET_SHAPES_DERIVATIVES_HELPER, \
			   INTEREST_LIST)

#undef AKANTU_SPECIALIZE_GET_SHAPE_DERIVATIVES_HELPER
#undef GET_SHAPES_DERIVATIVES
#undef INTEREST_LIST

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
inline const Array<Real> & FEEngineTemplate<I, S, kind>::getShapesDerivatives(const ElementType & type,
									      const GhostType & ghost_type,
									      __attribute__((unused)) UInt id) const {

  return GetShapesDerivativesHelper<kind>::call(shape_functions, type, ghost_type, id);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template<ElementKind kind>
struct GetQuadraturePointsHelper { };

#define GET_QUADS(type)                                                 \
  ret = &(integrator. template getQuadraturePoints<type>(ghost_type));	\

#define AKANTU_SPECIALIZE_GET_QUADRATURE_POINTS_HELPER(kind)		\
  template<>								\
  struct GetQuadraturePointsHelper<kind> {				\
    template <template <ElementKind> class I,				\
	      ElementKind k>						\
    static const Matrix<Real> & call(const I<k> & integrator,		\
				     const ElementType type,		\
				     const GhostType & ghost_type) {	\
      const Matrix<Real> * ret = NULL;					\
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(GET_QUADS, kind);		\
      return *ret;							\
    }									\
    };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_GET_QUADRATURE_POINTS_HELPER)

#undef AKANTU_SPECIALIZE_GET_QUADRATURE_POINTS_HELPER
#undef GET_QUADS

template<template <ElementKind> class I,
	 template <ElementKind> class S,
	 ElementKind kind>
inline const Matrix<Real> &
FEEngineTemplate<I, S, kind>::getQuadraturePoints(const ElementType & type,
						  const GhostType & ghost_type) const {
  return GetQuadraturePointsHelper<kind>::call(integrator, type, ghost_type);
}

/* -------------------------------------------------------------------------- */
__END_AKANTU__
#include "shape_lagrange.hh"
#include "integrator_gauss.hh"
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <>
template <>
inline void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular>::
assembleLumpedTemplate<_triangle_6>(const Array<Real> & field_1,
				    UInt nb_degree_of_freedom,
				    Array<Real> & lumped,
				    const Array<Int> & equation_number,
				    const GhostType & ghost_type) const {
  assembleLumpedDiagonalScaling<_triangle_6>(field_1, nb_degree_of_freedom,
					     lumped, equation_number, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
inline void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular>::
assembleLumpedTemplate<_tetrahedron_10>(const Array<Real> & field_1,
					UInt nb_degree_of_freedom,
					Array<Real> & lumped,
					const Array<Int> & equation_number,
					const GhostType & ghost_type) const {
  assembleLumpedDiagonalScaling<_tetrahedron_10>(field_1, nb_degree_of_freedom,
						 lumped,  equation_number, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
inline void
FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular>::
assembleLumpedTemplate<_quadrangle_8>(const Array<Real> & field_1,
				      UInt nb_degree_of_freedom,
				      Array<Real> & lumped,
				      const Array<Int> & equation_number,
				      const GhostType & ghost_type) const {
  assembleLumpedDiagonalScaling<_quadrangle_8>(field_1, nb_degree_of_freedom,
					       lumped, equation_number, ghost_type);
}

/* -------------------------------------------------------------------------- */
template<>
template<>
inline void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular>::
computeNormalsOnControlPoints<_point_1>(__attribute__((unused))const Array<Real> & field,
					Array<Real> & normal,
					const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(mesh.getSpatialDimension() == 1, "Mesh dimension must be 1 to compute normals on points!");
  const ElementType type = _point_1;
  UInt spatial_dimension = mesh.getSpatialDimension();
  //UInt nb_nodes_per_element  = Mesh::getNbNodesPerElement(type);
  UInt nb_points = getNbQuadraturePoints(type, ghost_type);

  UInt nb_element = mesh.getConnectivity(type, ghost_type).getSize();
  normal.resize(nb_element * nb_points);
  Array<Real>::matrix_iterator normals_on_quad = normal.begin_reinterpret(spatial_dimension,
									  nb_points,
									  nb_element);
  Array< std::vector<Element> > segments = mesh.getElementToSubelement(type, ghost_type);
  Array<Real> coords = mesh.getNodes();

  const Mesh * mesh_segment;
  if (mesh.isMeshFacets())
    mesh_segment = &(mesh.getMeshParent());
  else
    mesh_segment = &mesh;

  for (UInt elem = 0; elem < nb_element; ++elem) {
    UInt nb_segment = segments(elem).size();

    AKANTU_DEBUG_ASSERT(nb_segment > 0, "Impossible to compute a normal on a point connected to 0 segments");

    Real normal_value = 1;
    if (nb_segment == 1) {
      Element segment = segments(elem)[0];
      const Array<UInt> & segment_connectivity = mesh_segment->getConnectivity(segment.type, segment.ghost_type);
      //const Vector<UInt> & segment_points = segment_connectivity.begin(Mesh::getNbNodesPerElement(segment.type))[segment.element];
      Real difference;
      if(segment_connectivity(0) == elem) {
	difference = coords(elem)-coords(segment_connectivity(1));
      } else {
	difference = coords(elem)-coords(segment_connectivity(0));
      }
      normal_value = difference/std::abs(difference);
    }

    for(UInt n(0); n < nb_points; ++n) {
      (*normals_on_quad)(0, n) = normal_value;
    }
    ++normals_on_quad;
  }

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
