/**
 * @file   material_core_includes.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Wed Feb 03 2016
 *
 * @brief  List of materials for core package
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_CORE_INCLUDES_HH__
#define __AKANTU_MATERIAL_CORE_INCLUDES_HH__

/* -------------------------------------------------------------------------- */
/* Material list                                                              */
/* -------------------------------------------------------------------------- */
#ifndef AKANTU_CMAKE_LIST_MATERIALS

// elastic materials
#include "material_elastic.hh"
#include "material_elastic_linear_anisotropic.hh"
#include "material_elastic_orthotropic.hh"
#include "material_neohookean.hh"

// visco-elastic materials
#include "material_standard_linear_solid_deviatoric.hh"

// damage laws
#include "material_marigo.hh"
#include "material_mazars.hh"

// small-deformation plasticity
#include "material_linear_isotropic_hardening.hh"

#endif

#define AKANTU_CORE_MATERIAL_LIST                                              \
  ((2, (elastic, MaterialElastic)))((2, (neohookean, MaterialNeohookean)))(    \
      (2, (elastic_orthotropic, MaterialElasticOrthotropic)))(                 \
      (2, (elastic_anisotropic, MaterialElasticLinearAnisotropic)))(           \
      (2, (sls_deviatoric, MaterialStandardLinearSolidDeviatoric)))(           \
      (2, (marigo, MaterialMarigo)))((2, (mazars, MaterialMazars)))(           \
      (2, (plastic_linear_isotropic_hardening,                                 \
           MaterialLinearIsotropicHardening)))

#endif /* __AKANTU_MATERIAL_CORE_INCLUDES_HH__ */
