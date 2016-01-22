#===============================================================================
# @file   contact.cmake
#
# @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Mon Jan 18 2016
#
# @brief  package description for contact
#
# @section LICENSE
#
# Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
# Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
# Solides)
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

package_declare(contact
  DESCRIPTION "Use Contact package of Akantu"
  DEPENDS cpparray implicit optimization core_cxx11)


package_declare_sources(contact
  common/aka_point.hh
  common/aka_ball.cc
  common/aka_plane.hh
  common/aka_polytope.hh
  common/aka_ball.hh
  common/aka_timer.hh
  common/aka_tree.hh
  common/aka_bounding_box.hh
  common/aka_bounding_box.cc
  common/aka_geometry.hh
  common/aka_geometry.cc
  model/solid_mechanics/solid_mechanics_model_element.hh

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
  )


if(AKANTU_CONTACT)
  list(APPEND AKANTU_BOOST_COMPONENTS
    chrono
    system
    )
endif()

package_declare_documentation_files(contact
  manual-contact.tex
  figures/hertz_3D.png
  )

package_declare_documentation(contact
  "This package enables the contact mechanics engine for Akantu"
  )

package_set_package_system_dependency(contact deb-src
  libboost-system-dev libboost-chrono-dev)
package_set_package_system_dependency(contact deb
  libboost-system1.58 libboost-chrono1.58)
