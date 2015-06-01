%{
  #include "aka_common.hh"
%}

namespace akantu {
  %ignore getStaticParser;
  %ignore getUserParser;
  %ignore initialize(int & argc, char ** & argv);
  %ignore initialize(const std::string & input_file, int & argc, char ** & argv);
}

%typemap(in) (int argc, char *argv[]) {
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  $1 = PyList_Size($input);
  $2 = (char **) malloc(($1+1)*sizeof(char *));
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
   if ($2) free($2);
}

%inline %{
  namespace akantu {
    void initialize(const std::string & input_file) {
      int argc = 0;
      char ** argv = NULL;
      initialize(input_file, argc, argv);
    }
    void initialize() {
      int argc = 0;
      char ** argv = NULL;
      initialize(argc, argv);
    }

    void _initializeWithArgv(const std::string & input_file, int argc, char *argv[]) {
      initialize(input_file, argc, argv);
    }
    
  }
%}

%pythoncode %{
  import sys as _aka_sys
  def initializeWithArgv(input_file):
    _initializeWithArgv(input_file, _aka_sys.argv)
%}
%include "aka_common.hh"
