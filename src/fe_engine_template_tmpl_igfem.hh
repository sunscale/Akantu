/**
 * @file   shape_igfem_inline_impl.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
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
#include "shape_igfem.hh"
#include "integrator_gauss_igfem.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_FE_ENGINE_TEMPLATE_TMPL_IGFEM_HH__
#define __AKANTU_FE_ENGINE_TEMPLATE_TMPL_IGFEM_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* compatibility functions */
/* -------------------------------------------------------------------------- */

template <>
inline void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_igfem>::initShapeFunctions(const Array<Real> & nodes,
											    const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it  = mesh.firstType(element_dimension, ghost_type, _ek_igfem);
  Mesh::type_iterator end = mesh.lastType(element_dimension, ghost_type, _ek_igfem);
  for(; it != end; ++it) {
    ElementType type = *it;
    integrator.initIntegrator(nodes, type, ghost_type);

#define INIT(_type)							\
    do { 								\
      const Matrix<Real> & all_quads =					\
	integrator.template getQuadraturePoints<_type>(ghost_type);	\
      const Matrix<Real> & quads_1 =						\
	integrator.template getQuadraturePoints<ElementClassProperty<_type>::sub_element_type_1>(ghost_type); \
      const Matrix<Real> & quads_2 =						\
	integrator.template getQuadraturePoints<ElementClassProperty<_type>::sub_element_type_2>(ghost_type); \
      shape_functions.initShapeFunctions(nodes,				\
					 all_quads,			\
					 quads_1,			\
					 quads_2,			\
					 _type,				\
					 ghost_type);			\
    } while(0)

  AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(INIT);
#undef INIT

  }

  AKANTU_DEBUG_OUT();
}

#endif /* __AKANTU_FE_ENGINE_TEMPLATE_TMPL_IGFEM_HH__ */

__END_AKANTU__
