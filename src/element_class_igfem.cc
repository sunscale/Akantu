/**
 * @file   element_class_igfem.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Common part of igfem element_classes
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

using std::sqrt;

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<> ElementType ElementClass<_igfem_triangle_3>::p1_type        = _triangle_3;
template<> ElementType ElementClass<_igfem_triangle_3>::facet_type[]   = {_segment_2};
/* -------------------------------------------------------------------------- */
__END_AKANTU__
