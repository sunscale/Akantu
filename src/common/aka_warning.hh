// Intel warnings
#if defined(__INTEL_COMPILER)

#  if defined(AKANTU_WARNING_IGNORE_UNUSED_PARAMETER)
#  endif

// Clang Warnings
#elif defined (__clang__) // test clang to be sure that when we test for gnu it is only gnu
#  pragma clang diagnostic push
#  if defined(AKANTU_WARNING_IGNORE_UNUSED_PARAMETER)
#    pragma clang diagnostic ignored "-Wunused-parameter"
#  endif

// GCC warnings
#elif (defined(__GNUC__) || defined(__GNUG__))
#  define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#  if GCC_VERSION > 40600
#    pragma GCC diagnostic push
#  endif
#  if defined(AKANTU_WARNING_IGNORE_UNUSED_PARAMETER)
#    pragma GCC diagnostic ignored "-Wunused-parameter"
#  endif
#endif
