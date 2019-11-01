#===============================================================================
# @file   core.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Mon Jan 18 2016
#
# @brief  package description for core
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
package_declare(core NOT_OPTIONAL
  DESCRIPTION "core package for Akantu"
  FEATURES_PUBLIC cxx_strong_enums cxx_defaulted_functions
                  cxx_deleted_functions cxx_auto_type cxx_decltype_auto
  FEATURES_PRIVATE cxx_lambdas cxx_nullptr cxx_range_for
		   cxx_delegating_constructors
  DEPENDS INTERFACE akantu_iterators Boost)

package_declare_sources(core
  common/aka_array.cc
  common/aka_array.hh
  common/aka_array_tmpl.hh
  common/aka_array_printer.hh
  common/aka_bbox.hh
  common/aka_blas_lapack.hh
  common/aka_circular_array.hh
  common/aka_circular_array_inline_impl.cc
  common/aka_common.cc
  common/aka_common.hh
  common/aka_common_inline_impl.cc
  common/aka_csr.hh
  common/aka_element_classes_info_inline_impl.cc
  common/aka_enum_macros.hh
  common/aka_error.cc
  common/aka_error.hh
  common/aka_event_handler_manager.hh
  common/aka_extern.cc
  common/aka_factory.hh
  common/aka_fwd.hh
  common/aka_grid_dynamic.hh
  common/aka_math.cc
  common/aka_math.hh
  common/aka_math_tmpl.hh
  common/aka_memory.cc
  common/aka_memory.hh
  common/aka_memory_inline_impl.cc
  common/aka_named_argument.hh
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
  common/aka_voigthelper_tmpl.hh
  common/aka_voigthelper.cc
  common/aka_warning.hh
  common/aka_warning_restore.hh
  
  fe_engine/element_class.hh
  fe_engine/element_class_tmpl.hh
  fe_engine/element_classes/element_class_hexahedron_8_inline_impl.cc
  fe_engine/element_classes/element_class_hexahedron_20_inline_impl.cc
  fe_engine/element_classes/element_class_pentahedron_6_inline_impl.cc
  fe_engine/element_classes/element_class_pentahedron_15_inline_impl.cc
  fe_engine/element_classes/element_class_point_1_inline_impl.cc
  fe_engine/element_classes/element_class_quadrangle_4_inline_impl.cc
  fe_engine/element_classes/element_class_quadrangle_8_inline_impl.cc
  fe_engine/element_classes/element_class_segment_2_inline_impl.cc
  fe_engine/element_classes/element_class_segment_3_inline_impl.cc
  fe_engine/element_classes/element_class_tetrahedron_10_inline_impl.cc
  fe_engine/element_classes/element_class_tetrahedron_4_inline_impl.cc
  fe_engine/element_classes/element_class_triangle_3_inline_impl.cc
  fe_engine/element_classes/element_class_triangle_6_inline_impl.cc
  fe_engine/element_type_conversion.hh

  fe_engine/fe_engine.cc
  fe_engine/fe_engine.hh
  fe_engine/fe_engine_inline_impl.cc
  fe_engine/fe_engine_template.hh
  fe_engine/fe_engine_template_tmpl_field.hh
  fe_engine/fe_engine_template_tmpl.hh
  fe_engine/geometrical_element_property.hh
  fe_engine/geometrical_element_property.cc
  fe_engine/gauss_integration.cc
  fe_engine/gauss_integration_tmpl.hh
  fe_engine/integrator.hh
  fe_engine/integrator_gauss.hh
  fe_engine/integrator_gauss_inline_impl.cc
  fe_engine/interpolation_element_tmpl.hh
  fe_engine/integration_point.hh
  fe_engine/shape_functions.hh
  fe_engine/shape_functions.cc
  fe_engine/shape_functions_inline_impl.cc
  fe_engine/shape_lagrange_base.cc
  fe_engine/shape_lagrange_base.hh
  fe_engine/shape_lagrange_base_inline_impl.cc
  fe_engine/shape_lagrange.hh
  fe_engine/shape_lagrange_inline_impl.cc
  fe_engine/element.hh

  io/dumper/dumpable.hh
  io/dumper/dumpable.cc
  io/dumper/dumpable_dummy.hh
  io/dumper/dumpable_inline_impl.hh
  io/dumper/dumper_field.hh
  io/dumper/dumper_material_padders.hh
  io/dumper/dumper_filtered_connectivity.hh
  io/dumper/dumper_element_partition.hh

  io/mesh_io.cc
  io/mesh_io.hh
  io/mesh_io/mesh_io_diana.cc
  io/mesh_io/mesh_io_diana.hh
  io/mesh_io/mesh_io_msh.cc
  io/mesh_io/mesh_io_msh.hh
  #io/model_io.cc
  #io/model_io.hh

  io/parser/algebraic_parser.hh
  io/parser/input_file_parser.hh
  io/parser/parsable.cc
  io/parser/parsable.hh
  io/parser/parser.cc
  io/parser/parser_real.cc
  io/parser/parser_random.cc
  io/parser/parser_types.cc
  io/parser/parser_input_files.cc
  io/parser/parser.hh
  io/parser/parser_tmpl.hh
  io/parser/parser_grammar_tmpl.hh
  io/parser/cppargparse/cppargparse.hh
  io/parser/cppargparse/cppargparse.cc
  io/parser/cppargparse/cppargparse_tmpl.hh

  io/parser/parameter_registry.cc
  io/parser/parameter_registry.hh
  io/parser/parameter_registry_tmpl.hh

  mesh/element_group.cc
  mesh/element_group.hh
  mesh/element_group_inline_impl.cc
  mesh/element_type_map.cc
  mesh/element_type_map.hh
  mesh/element_type_map_tmpl.hh
  mesh/element_type_map_filter.hh
  mesh/group_manager.cc
  mesh/group_manager.hh
  mesh/group_manager_inline_impl.cc
  mesh/mesh.cc
  mesh/mesh.hh
  mesh/mesh_periodic.cc
  mesh/mesh_accessor.hh
  mesh/mesh_events.hh
  mesh/mesh_filter.hh
  mesh/mesh_global_data_updater.hh
  mesh/mesh_data.cc
  mesh/mesh_data.hh
  mesh/mesh_data_tmpl.hh
  mesh/mesh_inline_impl.cc
  mesh/node_group.cc
  mesh/node_group.hh
  mesh/node_group_inline_impl.cc
  mesh/mesh_iterators.hh

  mesh_utils/mesh_partition.cc
  mesh_utils/mesh_partition.hh
  mesh_utils/mesh_partition/mesh_partition_mesh_data.cc
  mesh_utils/mesh_partition/mesh_partition_mesh_data.hh
  mesh_utils/mesh_partition/mesh_partition_scotch.hh
  mesh_utils/mesh_utils_pbc.cc
  mesh_utils/mesh_utils.cc
  mesh_utils/mesh_utils.hh
  mesh_utils/mesh_utils_distribution.cc
  mesh_utils/mesh_utils_distribution.hh
  mesh_utils/mesh_utils.hh
  mesh_utils/mesh_utils_inline_impl.cc
  mesh_utils/global_ids_updater.hh
  mesh_utils/global_ids_updater.cc
  mesh_utils/global_ids_updater_inline_impl.cc

  model/common/boundary_condition/boundary_condition.hh
  model/common/boundary_condition/boundary_condition_functor.hh
  model/common/boundary_condition/boundary_condition_functor_inline_impl.cc
  model/common/boundary_condition/boundary_condition_tmpl.hh

  model/common/non_local_toolbox/neighborhood_base.hh
  model/common/non_local_toolbox/neighborhood_base.cc
  model/common/non_local_toolbox/neighborhood_base_inline_impl.cc
  model/common/non_local_toolbox/neighborhoods_criterion_evaluation/neighborhood_max_criterion.hh
  model/common/non_local_toolbox/neighborhoods_criterion_evaluation/neighborhood_max_criterion.cc
  model/common/non_local_toolbox/neighborhoods_criterion_evaluation/neighborhood_max_criterion_inline_impl.cc
  model/common/non_local_toolbox/non_local_manager.hh
  model/common/non_local_toolbox/non_local_manager.cc
  model/common/non_local_toolbox/non_local_manager_inline_impl.cc
  model/common/non_local_toolbox/non_local_manager_callback.hh
  model/common/non_local_toolbox/non_local_neighborhood_base.hh
  model/common/non_local_toolbox/non_local_neighborhood_base.cc
  model/common/non_local_toolbox/non_local_neighborhood.hh
  model/common/non_local_toolbox/non_local_neighborhood_tmpl.hh
  model/common/non_local_toolbox/non_local_neighborhood_inline_impl.cc
  model/common/non_local_toolbox/base_weight_function.hh
  model/common/non_local_toolbox/base_weight_function_inline_impl.cc

  model/common/model_solver.cc
  model/common/model_solver.hh
  model/common/solver_callback.hh
  model/common/solver_callback.cc

  model/common/dof_manager/dof_manager.cc
  model/common/dof_manager/dof_manager.hh
  model/common/dof_manager/dof_manager_default.cc
  model/common/dof_manager/dof_manager_default.hh
  model/common/dof_manager/dof_manager_default_inline_impl.cc
  model/common/dof_manager/dof_manager_inline_impl.cc

  model/common/non_linear_solver/non_linear_solver.cc
  model/common/non_linear_solver/non_linear_solver.hh
  model/common/non_linear_solver/non_linear_solver_default.hh
  model/common/non_linear_solver/non_linear_solver_lumped.cc
  model/common/non_linear_solver/non_linear_solver_lumped.hh

  model/common/time_step_solvers/time_step_solver.hh
  model/common/time_step_solvers/time_step_solver.cc
  model/common/time_step_solvers/time_step_solver_default.cc
  model/common/time_step_solvers/time_step_solver_default.hh
  model/common/time_step_solvers/time_step_solver_default_explicit.hh

  model/common/integration_scheme/generalized_trapezoidal.cc
  model/common/integration_scheme/generalized_trapezoidal.hh
  model/common/integration_scheme/integration_scheme.cc
  model/common/integration_scheme/integration_scheme.hh
  model/common/integration_scheme/integration_scheme_1st_order.cc
  model/common/integration_scheme/integration_scheme_1st_order.hh
  model/common/integration_scheme/integration_scheme_2nd_order.cc
  model/common/integration_scheme/integration_scheme_2nd_order.hh
  model/common/integration_scheme/newmark-beta.cc
  model/common/integration_scheme/newmark-beta.hh
  model/common/integration_scheme/pseudo_time.cc
  model/common/integration_scheme/pseudo_time.hh

  model/model.cc
  model/model.hh
  model/model_inline_impl.cc
  model/model_options.hh

  solver/solver_vector.hh
  solver/solver_vector_default.cc
  solver/solver_vector_default.hh
  solver/solver_vector_default_tmpl.hh
  solver/solver_vector_distributed.cc
  solver/solver_vector_distributed.hh
  solver/sparse_matrix.cc
  solver/sparse_matrix.hh
  solver/sparse_matrix_aij.cc
  solver/sparse_matrix_aij.hh
  solver/sparse_matrix_aij_inline_impl.cc
  solver/sparse_matrix_inline_impl.cc
  solver/sparse_solver.cc
  solver/sparse_solver.hh
  solver/sparse_solver_inline_impl.cc
  solver/terms_to_assemble.hh

  synchronizer/communication_buffer_inline_impl.cc
  synchronizer/communication_descriptor.hh
  synchronizer/communication_descriptor_tmpl.hh
  synchronizer/communication_request.hh
  synchronizer/communication_tag.hh
  synchronizer/communications.hh
  synchronizer/communications_tmpl.hh
  synchronizer/communicator.cc
  synchronizer/communicator.hh
  synchronizer/communicator_dummy_inline_impl.cc
  synchronizer/communicator_event_handler.hh
  synchronizer/communicator_inline_impl.hh
  synchronizer/data_accessor.cc
  synchronizer/data_accessor.hh
  synchronizer/dof_synchronizer.cc
  synchronizer/dof_synchronizer.hh
  synchronizer/dof_synchronizer_inline_impl.cc
  synchronizer/element_info_per_processor.cc
  synchronizer/element_info_per_processor.hh
  synchronizer/element_info_per_processor_tmpl.hh
  synchronizer/element_synchronizer.cc
  synchronizer/element_synchronizer.hh
  synchronizer/facet_synchronizer.cc
  synchronizer/facet_synchronizer.hh
  synchronizer/facet_synchronizer_inline_impl.cc
  synchronizer/grid_synchronizer.cc
  synchronizer/grid_synchronizer.hh
  synchronizer/grid_synchronizer_tmpl.hh
  synchronizer/master_element_info_per_processor.cc
  synchronizer/node_info_per_processor.cc
  synchronizer/node_info_per_processor.hh
  synchronizer/node_synchronizer.cc
  synchronizer/node_synchronizer.hh
  synchronizer/periodic_node_synchronizer.cc
  synchronizer/periodic_node_synchronizer.hh
  synchronizer/slave_element_info_per_processor.cc
  synchronizer/synchronizer.cc
  synchronizer/synchronizer.hh
  synchronizer/synchronizer_impl.hh
  synchronizer/synchronizer_impl_tmpl.hh
  synchronizer/synchronizer_registry.cc
  synchronizer/synchronizer_registry.hh
  synchronizer/synchronizer_tmpl.hh
  synchronizer/communication_buffer.hh
  )

set(AKANTU_SPIRIT_SOURCES
  io/mesh_io/mesh_io_abaqus.cc
  io/parser/parser_real.cc
  io/parser/parser_random.cc
  io/parser/parser_types.cc
  io/parser/parser_input_files.cc
  PARENT_SCOPE
  )

package_declare_elements(core
  ELEMENT_TYPES
  _point_1
  _segment_2
  _segment_3
  _triangle_3
  _triangle_6
  _quadrangle_4
  _quadrangle_8
  _tetrahedron_4
  _tetrahedron_10
  _pentahedron_6
  _pentahedron_15
  _hexahedron_8
  _hexahedron_20
  KIND regular
  GEOMETRICAL_TYPES
  _gt_point
  _gt_segment_2
  _gt_segment_3
  _gt_triangle_3
  _gt_triangle_6
  _gt_quadrangle_4
  _gt_quadrangle_8
  _gt_tetrahedron_4
  _gt_tetrahedron_10
  _gt_hexahedron_8
  _gt_hexahedron_20
  _gt_pentahedron_6
  _gt_pentahedron_15
  INTERPOLATION_TYPES
  _itp_lagrange_point_1
  _itp_lagrange_segment_2
  _itp_lagrange_segment_3
  _itp_lagrange_triangle_3
  _itp_lagrange_triangle_6
  _itp_lagrange_quadrangle_4
  _itp_serendip_quadrangle_8
  _itp_lagrange_tetrahedron_4
  _itp_lagrange_tetrahedron_10
  _itp_lagrange_hexahedron_8
  _itp_serendip_hexahedron_20
  _itp_lagrange_pentahedron_6
  _itp_lagrange_pentahedron_15
  GEOMETRICAL_SHAPES
  _gst_point
  _gst_triangle
  _gst_square
  _gst_prism
  GAUSS_INTEGRATION_TYPES
  _git_point
  _git_segment
  _git_triangle
  _git_tetrahedron
  _git_pentahedron
  INTERPOLATION_KIND _itk_lagrangian
  FE_ENGINE_LISTS
  gradient_on_integration_points
  interpolate_on_integration_points
  interpolate
  compute_normals_on_integration_points
  inverse_map
  contains
  compute_shapes
  compute_shapes_derivatives
  get_shapes_derivatives
  lagrange_base
  )

package_declare_documentation_files(core
  manual.sty
  manual.cls
  manual.tex
  manual-macros.sty
  manual-titlepages.tex
  manual-authors.tex
  manual-changelog.tex
  manual-introduction.tex
  manual-gettingstarted.tex
  manual-io.tex
  manual-feengine.tex
  manual-elements.tex
  manual-appendix-elements.tex
  manual-appendix-packages.tex
  manual-backmatter.tex
  manual-bibliography.bib
  manual-bibliographystyle.bst

  figures/bc_and_ic_example.pdf
  figures/boundary.pdf
  figures/boundary.svg
  figures/dirichlet.pdf
  figures/dirichlet.svg
#  figures/doc_wheel.pdf
#  figures/doc_wheel.svg
  figures/hot-point-1.png
  figures/hot-point-2.png
  figures/insertion.pdf
  figures/interpolate.pdf
  figures/interpolate.svg
  figures/vectors.pdf
  figures/vectors.svg

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

package_declare_documentation(core
  "This package is the core engine of \\akantu. It depends on:"
  "\\begin{itemize}"
  "\\item A C++ compiler (\\href{http://gcc.gnu.org/}{GCC} >= 4, or \\href{https://software.intel.com/en-us/intel-compilers}{Intel})."
  "\\item The cross-platform, open-source \\href{http://www.cmake.org/}{CMake} build system."
  "\\item The \\href{http://www.boost.org/}{Boost} C++ portable libraries."
  "\\item The \\href{http://www.zlib.net/}{zlib} compression library."
  "\\end{itemize}"
  ""
  "Under Ubuntu (14.04 LTS) the installation can be performed using the commands:"
  "\\begin{command}"
  "  > sudo apt-get install cmake libboost-dev zlib1g-dev g++"
  "\\end{command}"
  ""
  "Under Mac OS X the installation requires the following steps:"
  "\\begin{itemize}"
  "\\item Install Xcode"
  "\\item Install the command line tools."
  "\\item Install the MacPorts project which allows to automatically"
  "download and install opensource packages."
  "\\end{itemize}"
  "Then the following commands should be typed in a terminal:"
  "\\begin{command}"
  "  > sudo port install cmake gcc48 boost"
  "\\end{command}"
  )

find_program(READLINK_COMMAND readlink)
find_program(ADDR2LINE_COMMAND addr2line)
find_program(PATCH_COMMAND patch)
mark_as_advanced(READLINK_COMMAND)
mark_as_advanced(ADDR2LINE_COMMAND)

package_declare_extra_files_to_package(core
  SOURCES
    common/aka_element_classes_info.hh.in
    common/aka_config.hh.in
  )

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.9))
  package_set_compile_flags(core CXX "-Wno-undefined-var-template")
endif()

if(DEFINED AKANTU_CXX11_FLAGS)
  package_declare(core_cxx11 NOT_OPTIONAL
    DESCRIPTION "C++ 11 additions for Akantu core"
    COMPILE_FLAGS CXX "${AKANTU_CXX11_FLAGS}")

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.6")
      set(AKANTU_CORE_CXX11 OFF CACHE BOOL "C++ 11 additions for Akantu core - not supported by the selected compiler" FORCE)
    endif()
  endif()

  package_declare_documentation(core_cxx11
    "This option activates some features of the C++11 standard. This is usable with GCC>=4.7 or Intel>=13.")
else()
  if(CMAKE_VERSION VERSION_LESS 3.1)
    message(FATAL_ERROR "Since version 3.0 Akantu requires at least c++11 capable compiler")
  endif()
endif()
