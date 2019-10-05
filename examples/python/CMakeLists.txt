add_subdirectory(custom-material)
add_subdirectory(dynamics)
add_subdirectory(eigen_modes)
add_subdirectory(plate-hole)
add_subdirectory(stiffness_matrix)

package_add_files_to_package(
  examples/python/README.rst
  )
