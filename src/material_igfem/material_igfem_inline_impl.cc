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
inline void MaterialIGFEM::interpolateInternal(const Element new_el,
					       const Element old_el,
					       Vector<Real> & interpolated,
					       const Vector<Real> & internal,
					       const UInt nb_quads_new,
					       const UInt nb_quads_old) {
  /// @todo make this function generic!! Works right now only for
  /// elements with constant fields! A generic function would map the
  /// sub element quads into the physical domain, then interpolate the
  /// _ek_regular element on these points.
  
  for (UInt q = 0; q < nb_quads_new; ++q)
    interpolated(q) = internal(0);
}


