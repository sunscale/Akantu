/**
 * @file   material_igfem_inline_impl.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief Implementation of the inline functions of the parent
 * material for IGFEM
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

__END_AKANTU__

#include "solid_mechanics_model_igfem.hh"
#include "igfem_helper.hh"
#include <iostream>

__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
/// map internals from an old regular element to new IGFEM element
inline void MaterialIGFEM::interpolateInternal(const ElementType & type,
				   const Vector<Real> & internal,
				   Vector<Real> & interpolated,
				   const UInt nb_quads,
				   const UInt sub_element) {
  /// @todo make this function generic!! Works right now only for
  /// elements with constant fields! A generic function would map the
  /// sub element quads into the physical domain, then interpolate the
  /// _ek_regular element on these points. For this operation the
  /// element coordinates would be needed as a parameter of the
  /// function
  
  UInt nb_quads_sub_element = IGFEMHelper::getNbQuadraturePoints(type, sub_element);
  UInt components = internal.size() / nb_quads;
  UInt start_idx = 0;
  if (sub_element) 
    start_idx = (nb_quads - nb_quads_sub_element) * components;
  for (UInt i = 0; i < nb_quads_sub_element; ++i)
    interpolated(start_idx + i) = internal(0);
}


