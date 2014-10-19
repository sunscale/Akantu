#===============================================================================
# @file   00_core.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Fri Sep 19 2014
#
# @brief  package description for core
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

set(AKANTU_CORE ON CACHE INTERNAL "core package for Akantu" FORCE)

set(AKANTU_CORE_FILES
  common/aka_array.cc
  common/aka_array.hh
  common/aka_array_tmpl.hh
  common/aka_blas_lapack.hh
  common/aka_circular_array.hh
  common/aka_circular_array_inline_impl.cc
  common/aka_common.cc
  common/aka_common.hh
  common/aka_common_inline_impl.cc
  common/aka_csr.hh
  common/aka_error.cc
  common/aka_error.hh
  common/aka_event_handler_manager.hh
  common/aka_extern.cc
  common/aka_fwd.hh
  common/aka_grid_dynamic.hh
  common/aka_math.cc
  common/aka_math.hh
  common/aka_math_tmpl.hh
  common/aka_memory.cc
  common/aka_memory.hh
  common/aka_memory_inline_impl.cc
  common/aka_random_generator.hh
  common/aka_safe_enum.hh
  common/aka_static_memory.cc
  common/aka_static_memory.hh
  common/aka_static_memory_inline_impl.cc
  common/aka_static_memory_tmpl.hh
  common/aka_typelist.hh
  common/aka_types.hh
  common/aka_visitor.hh
  common/aka_voigthelper.hh
  common/aka_voigthelper.cc

  fe_engine/element_class.cc
  fe_engine/element_class.hh
  fe_engine/element_class_tmpl.hh
  fe_engine/element_classes/element_class_hexahedron_8_inline_impl.cc
  fe_engine/element_classes/element_class_pentahedron_6_inline_impl.cc
  fe_engine/element_classes/element_class_point_1_inline_impl.cc
  fe_engine/element_classes/element_class_quadrangle_4_inline_impl.cc
  fe_engine/element_classes/element_class_quadrangle_8_inline_impl.cc
  fe_engine/element_classes/element_class_segment_2_inline_impl.cc
  fe_engine/element_classes/element_class_segment_3_inline_impl.cc
  fe_engine/element_classes/element_class_tetrahedron_10_inline_impl.cc
  fe_engine/element_classes/element_class_tetrahedron_4_inline_impl.cc
  fe_engine/element_classes/element_class_triangle_3_inline_impl.cc
  fe_engine/element_classes/element_class_triangle_6_inline_impl.cc

  fe_engine/fe_engine.cc
  fe_engine/fe_engine.hh
  fe_engine/fe_engine_inline_impl.cc
  fe_engine/fe_engine_template.hh
  fe_engine/fe_engine_template_tmpl.hh
  fe_engine/geometrical_data_tmpl.hh
  fe_engine/geometrical_element.cc
  fe_engine/integration_element.cc
  fe_engine/integrator.hh
  fe_engine/integrator_gauss.hh
  fe_engine/integrator_gauss_inline_impl.cc
  fe_engine/interpolation_element.cc
  fe_engine/interpolation_element_tmpl.hh
  fe_engine/shape_functions.hh
  fe_engine/shape_functions_inline_impl.cc
  fe_engine/shape_lagrange.cc
  fe_engine/shape_lagrange.hh
  fe_engine/shape_lagrange_inline_impl.cc
  fe_engine/shape_linked.cc
  fe_engine/shape_linked.hh
  fe_engine/shape_linked_inline_impl.cc
  fe_engine/element.hh

  io/dumper/dumpable.hh
  io/dumper/dumpable.cc
  io/dumper/dumpable_inline_impl.hh
  io/dumper/dumper_field.hh
  io/dumper/dumper_material_padders.hh
  io/dumper/dumper_filtered_connectivity.hh
  io/dumper/dumper_element_type.hh
  io/dumper/dumper_element_partition.hh

  io/mesh_io.cc
  io/mesh_io.hh
  io/mesh_io/mesh_io_abaqus.cc
  io/mesh_io/mesh_io_abaqus.hh
  io/mesh_io/mesh_io_diana.cc
  io/mesh_io/mesh_io_diana.hh
  io/mesh_io/mesh_io_msh.cc
  io/mesh_io/mesh_io_msh.hh
  io/model_io.cc
  io/model_io.hh

  io/parser/algebraic_parser.hh
  io/parser/input_file_parser.hh
  io/parser/parsable.cc
  io/parser/parsable.hh
  io/parser/parsable_tmpl.hh
  io/parser/parser.cc
  io/parser/parser.hh
  io/parser/parser_tmpl.hh
  io/parser/cppargparse/cppargparse.hh
  io/parser/cppargparse/cppargparse.cc
  io/parser/cppargparse/cppargparse_tmpl.hh

  mesh/element_group.cc
  mesh/element_group.hh
  mesh/element_group_inline_impl.cc
  mesh/element_type_map.hh
  mesh/element_type_map_tmpl.hh
  mesh/element_type_map_filter.hh
  mesh/group_manager.cc
  mesh/group_manager.hh
  mesh/group_manager_inline_impl.cc
  mesh/mesh.cc
  mesh/mesh.hh
  mesh/mesh_filter.hh
  mesh/mesh_data.cc
  mesh/mesh_data.hh
  mesh/mesh_data_tmpl.hh
  mesh/mesh_inline_impl.cc
  mesh/node_group.cc
  mesh/node_group.hh
  mesh/node_group_inline_impl.cc

  mesh_utils/mesh_partition.cc
  mesh_utils/mesh_partition.hh
  mesh_utils/mesh_partition/mesh_partition_mesh_data.cc
  mesh_utils/mesh_partition/mesh_partition_mesh_data.hh
  mesh_utils/mesh_partition/mesh_partition_scotch.hh
  mesh_utils/mesh_pbc.cc
  mesh_utils/mesh_utils.cc
  mesh_utils/mesh_utils.hh
  mesh_utils/mesh_utils_inline_impl.cc

  model/boundary_condition.hh
  model/boundary_condition_functor.hh
  model/boundary_condition_functor_inline_impl.cc
  model/boundary_condition_tmpl.hh
  model/integration_scheme/generalized_trapezoidal.hh
  model/integration_scheme/generalized_trapezoidal_inline_impl.cc
  model/integration_scheme/integration_scheme_1st_order.hh
  model/integration_scheme/integration_scheme_2nd_order.hh
  model/integration_scheme/newmark-beta.hh
  model/integration_scheme/newmark-beta_inline_impl.cc
  model/model.cc
  model/model.hh
  model/model_inline_impl.cc

  model/solid_mechanics/material.cc
  model/solid_mechanics/material.hh
  model/solid_mechanics/material_inline_impl.cc
  model/solid_mechanics/material_list.hh
  model/solid_mechanics/material_random_internal.hh
  model/solid_mechanics/material_selector.hh
  model/solid_mechanics/material_selector_tmpl.hh
  model/solid_mechanics/materials/internal_field.hh
  model/solid_mechanics/materials/internal_field_tmpl.hh
  model/solid_mechanics/materials/random_internal_field.hh
  model/solid_mechanics/materials/random_internal_field_tmpl.hh
  model/solid_mechanics/solid_mechanics_model.cc
  model/solid_mechanics/solid_mechanics_model.hh
  model/solid_mechanics/solid_mechanics_model_inline_impl.cc
  model/solid_mechanics/solid_mechanics_model_mass.cc
  model/solid_mechanics/solid_mechanics_model_material.cc
  model/solid_mechanics/solid_mechanics_model_tmpl.hh
  model/solid_mechanics/solid_mechanics_model_event_handler.hh
  model/solid_mechanics/materials/plane_stress_toolbox.hh
  model/solid_mechanics/materials/plane_stress_toolbox_tmpl.hh


  model/solid_mechanics/materials/material_elastic.cc
  model/solid_mechanics/materials/material_elastic.hh
  model/solid_mechanics/materials/material_elastic_inline_impl.cc
  model/solid_mechanics/materials/material_thermal.cc
  model/solid_mechanics/materials/material_thermal.hh
  model/solid_mechanics/materials/material_elastic_linear_anisotropic.cc
  model/solid_mechanics/materials/material_elastic_linear_anisotropic.hh
  model/solid_mechanics/materials/material_elastic_orthotropic.cc
  model/solid_mechanics/materials/material_elastic_orthotropic.hh
  model/solid_mechanics/materials/material_damage/material_damage.hh
  model/solid_mechanics/materials/material_damage/material_damage_tmpl.hh
  model/solid_mechanics/materials/material_damage/material_marigo.cc
  model/solid_mechanics/materials/material_damage/material_marigo.hh
  model/solid_mechanics/materials/material_damage/material_marigo_inline_impl.cc
  model/solid_mechanics/materials/material_damage/material_mazars.cc
  model/solid_mechanics/materials/material_damage/material_mazars.hh
  model/solid_mechanics/materials/material_damage/material_mazars_inline_impl.cc
  model/solid_mechanics/materials/material_finite_deformation/material_neohookean.cc
  model/solid_mechanics/materials/material_finite_deformation/material_neohookean.hh
  model/solid_mechanics/materials/material_finite_deformation/material_neohookean_inline_impl.cc
  model/solid_mechanics/materials/material_plastic/material_plastic.cc
  model/solid_mechanics/materials/material_plastic/material_plastic.hh
  model/solid_mechanics/materials/material_plastic/material_plastic_inline_impl.cc
  model/solid_mechanics/materials/material_plastic/material_linear_isotropic_hardening.cc
  model/solid_mechanics/materials/material_plastic/material_linear_isotropic_hardening.hh
  model/solid_mechanics/materials/material_plastic/material_linear_isotropic_hardening_inline_impl.cc
  model/solid_mechanics/materials/material_viscoelastic/material_standard_linear_solid_deviatoric.cc
  model/solid_mechanics/materials/material_viscoelastic/material_standard_linear_solid_deviatoric.hh



  solver/solver.cc
  solver/solver.hh
  solver/sparse_matrix.cc
  solver/sparse_matrix.hh
  solver/sparse_matrix_inline_impl.cc
  solver/static_solver.hh
  solver/static_solver.cc

  synchronizer/communication_buffer.hh
  synchronizer/communication_buffer_inline_impl.cc
  synchronizer/data_accessor.cc
  synchronizer/data_accessor.hh
  synchronizer/data_accessor_inline_impl.cc
  synchronizer/distributed_synchronizer.cc
  synchronizer/distributed_synchronizer.hh
  synchronizer/distributed_synchronizer_tmpl.hh
  synchronizer/dof_synchronizer.cc
  synchronizer/dof_synchronizer.hh
  synchronizer/dof_synchronizer_inline_impl.cc
  synchronizer/filtered_synchronizer.cc
  synchronizer/filtered_synchronizer.hh
  synchronizer/mpi_type_wrapper.hh
  synchronizer/pbc_synchronizer.cc
  synchronizer/pbc_synchronizer.hh
  synchronizer/real_static_communicator.hh
  synchronizer/static_communicator.cc
  synchronizer/static_communicator.hh
  synchronizer/static_communicator_dummy.hh
  synchronizer/static_communicator_inline_impl.hh
  synchronizer/synchronizer.cc
  synchronizer/synchronizer.hh
  synchronizer/synchronizer_registry.cc
  synchronizer/synchronizer_registry.hh
  )


set(AKANTU_CORE_DEB_DEPEND
  libboost-dev
  )

set(AKANTU_CORE_TESTS
  test_csr
  test_facet_element_mapping
  test_facet_extraction_tetrahedron_4
  test_facet_extraction_triangle_3
  test_grid
  test_interpolate_stress
  test_local_material
  test_material_damage_non_local
  test_material_thermal
  test_matrix
  test_mesh_boundary
  test_mesh_data
  test_mesh_io_msh
  test_mesh_io_msh_physical_names
  test_mesh_partitionate_mesh_data
  test_parser
  test_dumper
  test_pbc_tweak
  test_purify_mesh
  test_solid_mechanics_model_bar_traction2d
  test_solid_mechanics_model_bar_traction2d_structured
  test_solid_mechanics_model_bar_traction2d_structured_pbc
  test_solid_mechanics_model_boundary_condition
  test_solid_mechanics_model_circle_2
  test_solid_mechanics_model_cube3d
  test_solid_mechanics_model_cube3d_pbc
  test_solid_mechanics_model_cube3d_tetra10
  test_solid_mechanics_model_square
  test_solid_mechanics_model_material_eigenstrain
  test_static_memory
  test_surface_extraction_tetrahedron_4
  test_surface_extraction_triangle_3
  test_vector
  test_vector_iterator
  test_weight
  test_math
  test_material_standard_linear_solid_deviatoric_relaxation
  test_material_standard_linear_solid_deviatoric_relaxation_tension
  test_material_plasticity
  )

set(AKANTU_CORE_MANUAL_FILES
  manual.sty
  manual.cls
  manual.tex
  manual-macros.sty
  manual-titlepages.tex
  manual-introduction.tex
  manual-gettingstarted.tex
  manual-io.tex
  manual-solidmechanicsmodel.tex
  manual-constitutive-laws.tex
  manual-lumping.tex
  manual-elements.tex
  manual-appendix-elements.tex
  manual-appendix-materials.tex
  manual-appendix-packages.tex
  manual-backmatter.tex
  manual-bibliography.bib
  manual-bibliographystyle.bst

  figures/bc_and_ic_example.pdf
  figures/boundary.pdf
  figures/boundary.svg
  figures/dirichlet.pdf
  figures/dirichlet.svg
  figures/doc_wheel.pdf
  figures/doc_wheel.svg
  figures/dynamic_analysis.png
  figures/explicit_dynamic.pdf
  figures/explicit_dynamic.svg
  figures/static.pdf
  figures/static.svg
  figures/hooke_law.pdf
  figures/hot-point-1.png
  figures/hot-point-2.png
  figures/implicit_dynamic.pdf
  figures/implicit_dynamic.svg
  figures/insertion.pdf
  figures/interpolate.pdf
  figures/interpolate.svg
  figures/problemDomain.pdf_tex
  figures/problemDomain.pdf
  figures/static_analysis.png
  figures/stress_strain_el.pdf
  figures/tangent.pdf
  figures/tangent.svg
  figures/vectors.pdf
  figures/vectors.svg

  figures/stress_strain_neo.pdf
  figures/visco_elastic_law.pdf
  figures/isotropic_hardening_plasticity.pdf
  figures/stress_strain_visco.pdf

  figures/elements/hexahedron_8.pdf
  figures/elements/hexahedron_8.svg
  figures/elements/quadrangle_4.pdf
  figures/elements/quadrangle_4.svg
  figures/elements/quadrangle_8.pdf
  figures/elements/quadrangle_8.svg
  figures/elements/segment_2.pdf
  figures/elements/segment_2.svg
  figures/elements/segment_3.pdf
  figures/elements/segment_3.svg
  figures/elements/tetrahedron_10.pdf
  figures/elements/tetrahedron_10.svg
  figures/elements/tetrahedron_4.pdf
  figures/elements/tetrahedron_4.svg
  figures/elements/triangle_3.pdf
  figures/elements/triangle_3.svg
  figures/elements/triangle_6.pdf
  figures/elements/triangle_6.svg
  figures/elements/xtemp.pdf
  )

find_program(READLINK_COMMAND readlink)
find_program(ADDR2LINE_COMMAND addr2line)
mark_as_advanced(READLINK_COMMAND)
mark_as_advanced(ADDR2LINE_COMMAND)

include(CheckFunctionExists)

check_function_exists(clock_gettime _clock_gettime)

if(NOT _clock_gettime)
  set(AKANTU_USE_OBSOLETE_GETTIMEOFDAY ON)
else()
  set(AKANTU_USE_OBSOLETE_GETTIMEOFDAY OFF)
endif()

set(AKANTU_CORE_DOCUMENTATION
"
This package is the core engine of \\akantu. It depends on:
\\begin{itemize}
\\item A C++ compiler (\\href{http://gcc.gnu.org/}{GCC} >= 4, or \\href{https://software.intel.com/en-us/intel-compilers}{Intel}).
\\item The cross-platform, open-source \\href{http://www.cmake.org/}{CMake} build system.
\\item The \\href{http://www.boost.org/}{Boost} C++ portable libraries.
\\item The \\href{http://www.zlib.net/}{zlib} compression library.
\\end{itemize}

Under Ubuntu (14.04 LTS) the installation can be performed using the commands:
\\begin{command}
  > sudo apt-get install build-essential cmake-curses-gui libboost-dev zlib1g-dev
\\end{command}

Under Mac OS X the installation requires the following steps:
\\begin{itemize}
\\item Install Xcode
\\item Install the command line tools.
\\item Install the MacPorts project which allows to automatically 
download and install opensource packages. 
\\end{itemize}
Then the following commands should be typed in a terminal:
\\begin{command}
  > sudo port install cmake gcc48 boost
\\end{command}
")