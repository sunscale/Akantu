#Profiling
set(CMAKE_CXX_FLAGS_PROFILING "-g -gdb3 -pg -DNDEBUG -DAKANTU_NDEBUG -O2"
  CACHE STRING "Flags used by the compiler during profiling builds")
set(CMAKE_C_FLAGS_PROFILING "-g -gdb3 -pg -DNDEBUG -DAKANTU_NDEBUG -O2"
  CACHE STRING "Flags used by the compiler during profiling builds")
set(CMAKE_Fortran_FLAGS_PROFILING "-g -gdb3 -pg -DNDEBUG -DAKANTU_NDEBUG -O2"
  CACHE STRING "Flags used by the compiler during profiling builds")
set(CMAKE_EXE_LINKER_FLAGS_PROFILING "-pg"
  CACHE STRING "Flags used by the linker during profiling builds")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILING "-pg"
  CACHE STRING "Flags used by the linker during profiling builds")


# Sanitize the code
if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.2") OR
    CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_SANITIZE "-g -gdb3 -O2 -fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer"
    CACHE STRING "Flags used by the compiler during sanitining builds")
  set(CMAKE_C_FLAGS_SANITIZE "-g -gdb3 -O2 -fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer"
    CACHE STRING "Flags used by the compiler during sanitining builds")
  set(CMAKE_Fortran_FLAGS_SANITIZE "-g -gdb3 -O2 -fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer"
    CACHE STRING "Flags used by the compiler during sanitining builds")
  set(CMAKE_EXE_LINKER_FLAGS_SANITIZE "-fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer"
    CACHE STRING "Flags used by the linker during sanitining builds")
  set(CMAKE_SHARED_LINKER_FLAGS_SANITIZE "-fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer"
    CACHE STRING "Flags used by the linker during sanitining builds")

  mark_as_advanced(CMAKE_CXX_FLAGS_PROFILING CMAKE_C_FLAGS_PROFILING
    CMAKE_Fortran_FLAGS_PROFILING CMAKE_EXE_LINKER_FLAGS_PROFILING
    CMAKE_SHARED_LINKER_FLAGS_PROFILING CMAKE_SHARED_LINKER_FLAGS_SANITIZE
    CMAKE_CXX_FLAGS_SANITIZE CMAKE_C_FLAGS_SANITIZE
    CMAKE_Fortran_FLAGS_SANITIZE CMAKE_EXE_LINKER_FLAGS_SANITIZE
    )
endif()
