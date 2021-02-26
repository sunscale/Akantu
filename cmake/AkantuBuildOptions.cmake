find_package(Git)
find_program(UUENCODE_EXECUTABLE uuencode)

function(git_version_info TARGET_NAME)
  cmake_parse_arguments(_gvi
    "ALL"
    "OUTPUT_FILE"
    "SOURCES"
    ${ARGN})

  if(NOT _gvi_SOURCES)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} ls-files --full-name
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      OUTPUT_STRIP_TRAILING_WHITESPACE
      OUTPUT_VARIABLE __git_files)

    string(REPLACE "\n" ";" __git_files "${__git_files}")
    foreach(_item IN LISTS __git_files)
      get_filename_component(_abs_item "${_item}" ABSOLUTE
        BASE_DIR "${PROJECT_SOURCE_DIR}")
      if(EXISTS "${_abs_item}" AND
          NOT IS_DIRECTORY "${_abs_item}" AND
          NOT IS_SYMLINK "${_abs_item}")
        list(APPEND _gvi_SOURCES "${_abs_item}")
      endif()
    endforeach()
  endif()
  if(_gvi_OUTPUT_FILE)
    set(OUTPUT_FILE "${_gvi_OUTPUT_FILE}")
  else()
    set(OUTPUT_FILE "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.hh")
  endif()

  configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/source_informations.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/source_informations.sh
    @ONLY)

  set(__stamp_file ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.stamp)
  add_custom_command(
    VERBATIM
    OUTPUT ${__stamp_file}
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/source_informations.sh
    COMMAND ${CMAKE_COMMAND} -E touch ${__stamp_file}
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    DEPENDS ${_gvi_SOURCES}
    COMMENT "Generating git version information for ${PROJECT_NAME}"
    )

  set_source_files_properties(${OUTPUT_FILE} PROPERTIES GENERATED TRUE)

  unset(_all)
  if(${_gvi_ALL})
    set(_all ALL)
  endif()

  add_custom_target(
    ${TARGET_NAME} ${_all}
    DEPENDS ${__stamp_file}
    SOURCES ${_gvi_SOURCES}
    )
endfunction()
