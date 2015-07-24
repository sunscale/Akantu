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

  %ignore CSR::begin;
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
%include "aka_element_classes_info.hh"
%include "element.hh"


%inline %{
namespace akantu{
  template <typename T>
  class CSRIterator{

  public:
  CSRIterator(CSR<T> & csr,UInt row) {
    this->it =  csr.begin(row);
    this->end =  csr.end(row);
  };
  
  ~CSRIterator(){
  };

  T & __next_cpp(){
    if (this->it == this->end) AKANTU_SILENT_EXCEPTION("StopIteration");
    T & ref = *(this->it);
    ++this->it;
    return ref;
  }
  
  private:
  
  typename CSR<T>::iterator it;
  typename CSR<T>::iterator end;
  };
}
%}

%extend akantu::CSRIterator<akantu::Element>
{
%insert("python") %{
    def __iter__(self):
       return self

    def __next__(self):
       try:
         return self.__next_cpp()
       except Exception as e:
         raise StopIteration


    def next(self):
       return self.__next__()

%}
}

%extend akantu::CSR<akantu::Element>
{
  akantu::CSRIterator<akantu::Element> row(akantu::UInt row){
    return akantu::CSRIterator<akantu::Element>(*$self,row);
  }
}

%include "aka_csr.hh"
namespace akantu {
  %template (CSRUInt) CSR<UInt>;
  %template (CSRElement) CSR<Element>;
  %template (CSRIteratorElement) CSRIterator<Element>;
 }

