/**
 * @file   material_igfem_includes.hh
 *
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  list of IGFEM materials
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_CMAKE_LIST_MATERIALS
#include "material_igfem.hh"
#include "material_igfem_elastic.hh"
#include "material_igfem_iterative_stiffness_reduction.hh"
#include "material_igfem_saw_tooth_damage.hh"
#endif

#define AKANTU_IGFEM_MATERIAL_LIST                                             \
  ((2, (igfem_elastic, MaterialIGFEMElastic)))(                                \
      (2, (igfem_saw_tooth_damage, MaterialIGFEMSawToothDamage)))(             \
      (2, (igfem_iterative_stiffness_reduction,                                \
           MaterialIGFEMIterativeStiffnessReduction)))
