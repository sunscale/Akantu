if(__CMAKE_DEBUG_MESSAGES)
  return()
endif()
set(__CMAKE_DEBUG_MESSAGES TRUE)

macro(cmake_register_debug_message_module module)
  set(_CMAKE_DEBUG_MESSAGE_MODULES ${CMAKE_DEBUG_MESSAGE_MODULES})
  list(APPEND _CMAKE_DEBUG_MESSAGE_MODULES ${module})
  set(CMAKE_DEBUG_MESSAGE_MODULES "${_CMAKE_DEBUG_MESSAGE_MODULES}"
    CACHE INTERNAL "List of modules handled by the debug messages system" FORCE)
endmacro()

macro(cmake_activate_debug_message)
  set(_default FALSE)
  if(ARGC EQUAL 0)
    set(_default TRUE)
  endif()

  foreach(_module ${CMAKE_DEBUG_MESSAGE_MODULES})
    set(CMAKE_DEBUG_MESSAGE_${_module} ${_default} CACHE INTERNAL "" FORCE)
  endforeach()

  foreach(_module ${ARGN})
    set(CMAKE_DEBUG_MESSAGE_${_module} TRUE CACHE INTERNAL "" FORCE)
  endforeach()
endmacro()


macro(cmake_deactivate_debug_message)
  foreach(_module ${CMAKE_DEBUG_MESSAGE_MODULES})
    if(CMAKE_DEBUG_MESSAGE_${_module})
      set(CMAKE_DEBUG_MESSAGE_${_module} FALSE CACHE INTERNAL "" FORCE)
    endif()
  endforeach()
endmacro()

macro(cmake_debug_message module)
  if(CMAKE_DEBUG_MESSAGE_${module})
    message("${PROJECT_NAME} - ${module}: ${ARGN}")
  endif()
endmacro()