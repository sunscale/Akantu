%{
  #include "aka_common.hh"
  #include "aka_csr.hh"
  #include "element.hh"
%}

namespace akantu {
  %ignore getStaticParser;
  %ignore getUserParser;
  %ignore initialize(int & argc, char ** & argv);
  %ignore initialize(const std::string & input_file, int & argc, char ** & argv);
  extern const Array<UInt> empty_filter;

}

%typemap(in) (int argc, char *argv[]) {
  int i = 0;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }

  $1 = PyList_Size($input);
  $2 = new char *[$1+1];

  for (i = 0; i < $1; i++) {
    PyObject *s = PyList_GetItem($input,i);
    if (!PyString_Check(s)) {
        free($2);
        PyErr_SetString(PyExc_ValueError, "List items must be strings");
        return NULL;
    }
    $2[i] = PyString_AsString(s);
  }
  $2[i] = 0;
}

%typemap(freearg) (int argc, char *argv[]) {
%#if defined(__INTEL_COMPILER)
//#pragma warning ( disable : 383 )
%#elif defined (__clang__) // test clang to be sure that when we test for gnu it is only gnu
%#elif (defined(__GNUC__) || defined(__GNUG__))
%#  if __cplusplus > 199711L
%#    pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
%#  endif
%#endif

  delete [] $2;

%#if defined(__INTEL_COMPILER)
//#pragma warning ( disable : 383 )
%#elif defined (__clang__) // test clang to be sure that when we test for gnu it is only gnu
%#elif (defined(__GNUC__) || defined(__GNUG__))
%#  if __cplusplus > 199711L
%#    pragma GCC diagnostic pop
%#  endif
%#endif
}

%inline %{
  namespace akantu {
#if defined(AKANTU_USE_MPI)
    const int MPI=1;
#else
    const int MPI=0;
#endif
    void _initializeWithArgv(const std::string & input_file, int argc, char *argv[]) {
      initialize(input_file, argc, argv);
    }
  }
%}

%pythoncode %{
  import sys as _aka_sys
  def initialize(input_file="", argv=_aka_sys.argv):
      if sys_as.modules[__name__].MPI == 1:
         print "Loading mpi4py"
         try:
           import mpi4py
         except ImportError:
           pass

      _initializeWithArgv(input_file, argv)

%}

%include "aka_config.hh"
%include "aka_common.hh"
%include "aka_element_classes_info.hh"
%include "element.hh"



