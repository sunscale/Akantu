#===============================================================================
# @file   10_structural_mechanics.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Mon Jul 07 2014
#
# @brief  package description for structural mechanics
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

option(AKANTU_STRUCTURAL_MECHANICS "Use Structural mechanics model package of Akantu" OFF)

add_internal_package_dependencies(STRUCTURAL_MECHANICS IMPLICIT)

set(AKANTU_STRUCTURAL_MECHANICS_FILES
  fe_engine/element_class_structural.hh
  fe_engine/element_classes/element_class_bernoulli_beam_inline_impl.cc
  fe_engine/fe_engine_template_tmpl_struct.hh
  fe_engine/element_classes/element_class_kirchhoff_shell_inline_impl.cc
  io/mesh_io/mesh_io_msh_struct.cc
  io/mesh_io/mesh_io_msh_struct.hh
  io/model_io/model_io_ibarras.cc
  io/model_io/model_io_ibarras.hh
  model/structural_mechanics/structural_mechanics_model.cc
  model/structural_mechanics/structural_mechanics_model.hh
  model/structural_mechanics/structural_mechanics_model_boundary.cc
  model/structural_mechanics/structural_mechanics_model_inline_impl.cc
  model/structural_mechanics/structural_mechanics_model_mass.cc
  )

set(AKANTU_STRUCTURAL_MECHANICS_DOC
  manual/manual-structuralmechanicsmodel.tex
  )

set(AKANTU_STRUCTURAL_MECHANICS_TESTS
  test_structural_mechanics_model_bernoulli_beam_2
  test_structural_mechanics_model_boundary_bernoulli_beam_2
  test_structural_mechanics_model_bernoulli_beam_2_exemple_1_1
  test_structural_mechanics_model_bernoulli_beam_2_complicated
  test_structural_mechanics_model_bernoulli_beam_2_exemple_1_1_y
  test_structural_mechanics_model_bernoulli_beam_3_exemple_1_1_xy
  test_structural_mechanics_model_bernoulli_beam_3_exemple_1_1_zy
  test_structural_mechanics_model_bernoulli_beam_3_local_force
  test_structural_mechanics_model_bernoulli_beam_3_exercice_12_10_13 
  test_structural_mechanics_model_kirchhoff_shell_patch_test_4_5_5
  test_structural_mechanics_model_bernoulli_beam_dynamics
  )

set(AKANTU_STRUCTURAL_MECHANICS_MANUAL_FILES
  manual-structuralmechanicsmodel.tex
  manual-structuralmechanicsmodel-elements.tex

  figures/beam_example.pdf
  figures/elements/bernoulli_2.pdf
  figures/elements/bernoulli_2.svg
  )
  
set(AKANTU_STRUCTURAL_MECHANICS_DOCUMENTATION "
This package activates the compilation for the Structural Mechanics engine of Akantu
")
