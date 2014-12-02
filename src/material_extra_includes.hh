/**
 * @file   material_extra_includes.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Extra list of materials
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_CMAKE_LIST_MATERIALS

// visco-elastic materials
#include "material_stiffness_proportional.hh"

// damage materials
#include "material_brittle.hh"
#include "material_damage_iterative.hh"
#include "material_damage_linear.hh"
#include "material_vreepeerlings.hh"

// plasticity
#include "material_viscoplastic.hh"

#endif

#define  AKANTU_EXTRA_MATERIAL_LIST                                     \
  ((2, (damage_linear      , MaterialDamageLinear                 )))   \
  ((2, (brittle            , MaterialBrittle                      )))   \
  ((2, (damage_iterative   , MaterialDamageIterative              )))   \
  ((2, (vreepeerlings      , MaterialVreePeerlings                )))   \
  ((2, (ve_stiffness_prop  , MaterialStiffnessProportional        )))   \
  ((2, (visco_plastic      , MaterialViscoPlastic                 )))
