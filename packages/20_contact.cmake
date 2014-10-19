#===============================================================================
# @file   20_contact.cmake
#
# @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Mon Sep 15 2014
#
# @brief  package description for contact
#
# @section LICENSE
#
# Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

option(AKANTU_CONTACT "Use Contact package of Akantu" OFF)

set(AKANTU_CONTACT_FILES
  #cc files
  contact/discretization.cc
  contact/element.cc
  contact/friction.cc
  contact/resolution/resolution_augmented_lagrangian.cc
  contact/scheme.cc
  contact/search.cc
  contact/surface.cc
  contact/zone.cc
  contact/contact_impl.cc
  model/model_manager.cc

  # include files
  contact/contact_common.hh
  contact/contact_manager.hh
  contact/discretization.hh
  contact/element.hh
  contact/friction.hh
  contact/resolution.hh
  contact/resolution/resolution_augmented_lagrangian.hh
  contact/scheme.hh
  contact/search.hh
  contact/surface.hh
  contact/zone.hh
  contact/contact_impl.hh
  model/model_manager.hh
  model/contact_manager0.hh
  )

add_external_package_dependencies(contact cpparray)
add_internal_package_dependencies(contact implicit)
add_internal_package_dependencies(contact optimization)

if(AKANTU_CONTACT)
  list(APPEND AKANTU_BOOST_COMPONENTS
    chrono
    system
    )
endif()

mark_as_advanced(AKANTU_USE_CPPARRAY)

set(AKANTU_CONTACT_TESTS
  test_hertz_2D
  test_hertz_3D
  test_offset_1slave
  test_offset_2slaves
  test_acurnier_2D_1
  test_acurnier_2D_2
  test_acurnier_2D_3
  test_acurnier_3D_1
  test_acurnier_3D_2
  test_acurnier_3D_3
  )

set(AKANTU_CONTACT_MANUAL_FILES
  manual-contact.tex

  figures/hertz_3D.png
)

set(AKANTU_CONTACT_DOCUMENTATION
"
This package enables the contact mechanics engine for Akantu
")