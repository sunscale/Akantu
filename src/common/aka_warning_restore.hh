#if defined(__INTEL_COMPILER)
//#pragma warning ( disable : 383 )
#elif defined(__clang__) // test clang to be sure that when we test for gnu it
                         // is only gnu
#pragma clang diagnostic pop
#elif defined(__GNUG__)
#if GCC_VERSION > 40600
#pragma GCC diagnostic pop
#else
#if defined(AKANTU_WARNING_IGNORE_UNUSED_PARAMETER)
#pragma GCC diagnostic warning "-Wunused-parameter"
#endif
#endif
#endif

#undef AKANTU_WARNING_IGNORE_UNUSED_PARAMETER
