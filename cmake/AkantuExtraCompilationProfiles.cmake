#Profiling
set(CMAKE_CXX_FLAGS_PROFILING "-g -ggdb3 -pg -DNDEBUG -DAKANTU_NDEBUG -O2"
  CACHE STRING "Flags used by the compiler during profiling builds")
set(CMAKE_C_FLAGS_PROFILING "-g -ggdb3 -pg -DNDEBUG -DAKANTU_NDEBUG -O2"
  CACHE STRING "Flags used by the compiler during profiling builds")
set(CMAKE_Fortran_FLAGS_PROFILING "-g -ggdb3 -pg -DNDEBUG -DAKANTU_NDEBUG -O2"
  CACHE STRING "Flags used by the compiler during profiling builds")
set(CMAKE_EXE_LINKER_FLAGS_PROFILING "-pg"
  CACHE STRING "Flags used by the linker during profiling builds")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILING "-pg"
  CACHE STRING "Flags used by the linker during profiling builds")

mark_as_advanced(CMAKE_CXX_FLAGS_PROFILING CMAKE_C_FLAGS_PROFILING
    CMAKE_Fortran_FLAGS_PROFILING CMAKE_EXE_LINKER_FLAGS_PROFILING
    CMAKE_SHARED_LINKER_FLAGS_PROFILING)
# Sanitize the code
if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.2") OR
    CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

  set(_sanitize "-g -ggdb3 -O2 -fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer -fsanitize-blacklist=${PROJECT_SOURCE_DIR}/cmake/sanitize-blacklist.txt")

  set(CMAKE_CXX_FLAGS_SANITIZE ${_sanitize}
    CACHE STRING "Flags used by the compiler during sanitizing builds")
  set(CMAKE_C_FLAGS_SANITIZE ${_sanitize}
    CACHE STRING "Flags used by the compiler during sanitizing builds")
  set(CMAKE_Fortran_FLAGS_SANITIZE ${_sanitize}
    CACHE STRING "Flags used by the compiler during sanitizing builds")
  set(CMAKE_EXE_LINKER_FLAGS_SANITIZE ${_sanitize}
    CACHE STRING "Flags used by the linker during sanitizing builds")
  set(CMAKE_SHARED_LINKER_FLAGS_SANITIZE ${_sanitize}
    CACHE STRING "Flags used by the linker during sanitizing builds")

  mark_as_advanced(CMAKE_SHARED_LINKER_FLAGS_SANITIZE
    CMAKE_CXX_FLAGS_SANITIZE CMAKE_C_FLAGS_SANITIZE
    CMAKE_Fortran_FLAGS_SANITIZE CMAKE_EXE_LINKER_FLAGS_SANITIZE
    )
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

  set(_sanitize "-g -ggdb3 -O2 -fPIE -fsanitize=memory -fsanitize-memory-track-origins -fsanitize-recover=all -fno-omit-frame-pointer -fsanitize-blacklist=${PROJECT_SOURCE_DIR}/cmake/sanitize-blacklist.txt")

  set(CMAKE_CXX_FLAGS_SANITIZEMEMORY ${_sanitize}
    CACHE STRING "Flags used by the compiler during sanitizing builds")
  set(CMAKE_C_FLAGS_SANITIZEMEMORY ${_sanitize}
    CACHE STRING "Flags used by the compiler during sanitizing builds")
  set(CMAKE_Fortran_FLAGS_SANITIZEMEMORY ${_sanitize}
    CACHE STRING "Flags used by the compiler during sanitizing builds")
  set(CMAKE_EXE_LINKER_FLAGS_SANITIZEMEMORY ${_sanitize}
    CACHE STRING "Flags used by the linker during sanitizing builds")
  set(CMAKE_SHARED_LINKER_FLAGS_SANITIZEMEMORY ${_sanitize}
    CACHE STRING "Flags used by the linker during sanitizing builds")

  mark_as_advanced(CMAKE_SHARED_LINKER_FLAGS_SANITIZEMEMORY
    CMAKE_CXX_FLAGS_SANITIZEMEMORY CMAKE_C_FLAGS_SANITIZEMEMORY
    CMAKE_Fortran_FLAGS_SANITIZEMEMORY CMAKE_EXE_LINKER_FLAGS_SANITIZEMEMORY
    )
endif()
