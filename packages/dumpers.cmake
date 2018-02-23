package_declare(dumpers
  DEFAULT ON
  DESCRIPTION "Dumpers for Akantu"
  DEPENDS INTERFACE IOHelper)

package_declare_sources(dumpers
  io/dumper/dumpable_iohelper.hh
  io/dumper/dumper_compute.hh
  io/dumper/dumper_element_iterator.hh
  io/dumper/dumper_elemental_field.hh
  io/dumper/dumper_generic_elemental_field.hh
  io/dumper/dumper_generic_elemental_field_tmpl.hh
  io/dumper/dumper_homogenizing_field.hh
  io/dumper/dumper_internal_material_field.hh
  io/dumper/dumper_iohelper.cc
  io/dumper/dumper_iohelper.hh
  io/dumper/dumper_iohelper_paraview.cc
  io/dumper/dumper_iohelper_paraview.hh
  io/dumper/dumper_nodal_field.hh
  io/dumper/dumper_padding_helper.hh
  io/dumper/dumper_quadrature_point_iterator.hh
  io/dumper/dumper_text.cc
  io/dumper/dumper_text.hh
  io/dumper/dumper_type_traits.hh
  io/dumper/dumper_variable.hh
  )

package_declare_documentation(dumpers
  "This package activates the IOHelper facilities withing Akantu. This is mandatory if you want to be able to output Paraview files"
  "as well as any Dumper within Akantu."
  )

package_declare_extra_files_to_package(dumpers
  PROJECT
    third-party/cmake/iohelper.cmake
    cmake/Modules/FindIOHelper.cmake
  )
