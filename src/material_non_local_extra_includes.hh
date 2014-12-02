/**
 * @file   material_non_local_extra_includes.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Non local materials includes
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_CMAKE_LIST_MATERIALS
#  include "material_vreepeerlings_non_local.hh"
#  include "material_brittle_non_local.hh"
#  include "material_damage_iterative_non_local.hh"
#endif

#define AKANTU_MATERIAL_WEIGHT_FUNCTION_TMPL_LIST			\
  ((stress_wf, (StressBasedWeightFunction  )))				\
  ((damage_wf, (DamagedWeightFunction      )))				\
  ((remove_wf, (RemoveDamagedWeightFunction)))				\
  ((base_wf,   (BaseWeightFunction         )))

#define AKANTU_POSSIBLE_DAMAGE_PARENT_MATERIALS                         \
  ((sls_deviatoric      , (RemoveDamagedWeightFunction)(MaterialStandardLinearSolidDeviatoric))) \
  ((neohookean_base_wf  , (BaseWeightFunction)(MaterialNeohookean                            ))) \
  ((neohookean_remove_wf, (RemoveDamagedWeightFunction)(MaterialNeohookean                   ))) \
  ((elastic             , (RemoveDamagedWeightFunction)(MaterialElastic                      )))

#define AKANTU_MATERIAL_VREEPEERLINGS_WEIGHT_FUNCTION_TMPL_LIST		\
  AKANTU_MATERIAL_WEIGHT_FUNCTION_TMPL_LIST				\
  ((removed_damrate_wf, RemoveDamagedWithDamageRateWeightFunction))

#define AKANTU_DAMAGE_NON_LOCAL_MATERIAL_EXTRA_LIST			\
  ((3, (vreepeerlings_non_local, MaterialVreePeerlingsNonLocal,		\
	AKANTU_POSSIBLE_DAMAGE_PARENT_MATERIALS)))                      \
  ((3, (brittle_non_local       , MaterialBrittleNonLocal,		\
	AKANTU_MATERIAL_WEIGHT_FUNCTION_TMPL_LIST)))                    \
  ((3, (damage_iterative_non_local       , MaterialDamageIterativeNonLocal, \
	AKANTU_MATERIAL_WEIGHT_FUNCTION_TMPL_LIST)))
