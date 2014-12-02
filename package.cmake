#===============================================================================
# @file   package.cmake
#
# @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
#
# @brief  package description for interface-enriched generalized IGFEM
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
#===============================================================================

package_declare(IGFEM
  DESCRIPTION "Use Interface-enriched generalized FEM")

package_declare_sources(igfem
  element_class_igfem.hh
  shape_igfem.hh
  shape_igfem_inline_impl.cc
  element_class_igfem.cc
  element_class_igfem_triangle_3_inline_impl.cc
  igfem_element.hh
  igfem_element.cc
  )

set(AKANTU_IGFEM_TESTS
  test_igfem_integrate
)

