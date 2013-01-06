# (ADD_PCH_RULE  _header_filename _src_list)
# Version 7/26/2010 4:55pm
#
# use this macro before "add_executable"
#
# _header_filename
#	header to make a .gch
#
# _src_list
#   the variable name (do not use ${..}) which contains a
#     a list of sources (a.cpp b.cpp c.cpp ...)
#  This macro will append a header file to it, then this src_list can be used in
#	"add_executable..."
#
#
# Now a .gch file should be generated and gcc should use it.
#		(add -Winvalid-pch to the cpp flags to verify)
#
# make clean should delete the pch file
#
# example : ADD_PCH_RULE(headers.h myprog_SRCS)


function(ADD_PCH_RULE _header_filename _src_list)
  set(_gch_filename "${CMAKE_CURRENT_BINARY_DIR}/${_header_filename}.gch")

  list(APPEND ${_src_list} ${_gch_filename})

  set(_args ${CMAKE_CXX_FLAGS} -D__aka_inline__=inline)

  get_filename_component(_gch_filename_path ${_gch_filename} PATH)
  file(MAKE_DIRECTORY ${_gch_filename_path})

  #  list(APPEND _args -c ${CMAKE_CURRENT_SOURCE_DIR}/${_header_filename})
  list(APPEND _args -c ${CMAKE_CURRENT_SOURCE_DIR}/${_header_filename} -o ${_gch_filename} -Winvalid-pch)

  get_directory_property(DIRINC INCLUDE_DIRECTORIES)
  foreach (_inc ${DIRINC})
    list(APPEND _args "-I" ${_inc})
  endforeach(_inc ${DIRINC})
  separate_arguments(_args)

  set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${_header_filename} PROPERTIES GENERATED 1)
  set_source_files_properties(${_gch_filename} PROPERTIES GENERATED 1)
  add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/${_header_filename}
    COMMAND ${CMAKE_COMMAND} -E copy  ${CMAKE_CURRENT_SOURCE_DIR}/${_header_filename} ${CMAKE_CURRENT_BINARY_DIR}/${_header_filename} # ensure same directory! Required by gcc
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_header_filename}
    )

  add_custom_command(OUTPUT ${_gch_filename}
    COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1} ${_args}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_header_filename}
    )
endfunction(ADD_PCH_RULE _header_filename _src_list)
