%{
#define SWIG_FILE_WITH_INIT
#include "aka_array.hh"
%}

%include "typemaps.i"

namespace akantu {
  %ignore Array::operator=;
  %ignore Array::operator[];
  %ignore Array::operator();
  %ignore Array::set;
  %ignore Array::begin;
  %ignore Array::end;
  %ignore Array::begin_reinterpret;
  %ignore Array::end_reinterpret;
};

%include "aka_array.hh"

namespace akantu {
  %template(RArray) Array<akantu::Real, true>;
  %template(UArray) Array<akantu::UInt, true>;
  %template(BArray) Array<bool, true>;
}



%include "numpy.i"
%init %{
  import_array();
%}

%define %akantu_array_typemaps(DATA_TYPE, DATA_TYPECODE)

%typemap(out, fragment="NumPy_Fragments") (akantu::Array< DATA_TYPE > &)
{
  npy_intp dims[2] = {$1->getSize(), $1->getNbComponent()};
  PyObject* obj = PyArray_SimpleNewFromData(2, dims, DATA_TYPECODE, $1->storage());
  PyArrayObject* array = (PyArrayObject*) obj;

  if (!array) SWIG_fail;
  $result = SWIG_Python_AppendOutput($result, obj);
}

%enddef


%akantu_array_typemaps(double,       NPY_DOUBLE)
%akantu_array_typemaps(float,        NPY_FLOAT )
%akantu_array_typemaps(unsigned int, NPY_UINT  )


%typemap(out, fragment="NumPy_Fragments") (akantu::Array< bool > &)
{
  npy_intp dims[2] = {$1->getSize(), $1->getNbComponent()};

  int data_typecode = NPY_NOTYPE;
  size_t s = sizeof(bool);
  switch(s) {
  case 1: data_typecode = NPY_BOOL; break;
  case 2: data_typecode = NPY_UINT8; break;
  case 4: data_typecode = NPY_UINT16; break;
  case 8: data_typecode = NPY_UINT32; break;
  }

  PyObject* obj = PyArray_SimpleNewFromData(2, dims, data_typecode, $1->storage());
  PyArrayObject* array = (PyArrayObject*) obj;

  if (!array) SWIG_fail;
  $result = SWIG_Python_AppendOutput($result, obj);
}

