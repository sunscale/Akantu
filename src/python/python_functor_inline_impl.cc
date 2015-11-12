#ifndef __AKANTU_PYTHON_FUNCTOR_INLINE_IMPL_CC__
#define __AKANTU_PYTHON_FUNCTOR_INLINE_IMPL_CC__
/* -------------------------------------------------------------------------- */
#include <numpy/arrayobject.h>
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
/* -------------------------------------------------------------------------- */



template <typename T>
inline int PythonFunctor::getPythonDataTypeCode() const{
  AKANTU_EXCEPTION("undefined type");
}

/* -------------------------------------------------------------------------- */

template <>
inline int PythonFunctor::getPythonDataTypeCode<bool>() const{
  int data_typecode = NPY_NOTYPE;
  size_t s = sizeof(bool);
  switch(s) {
  case 1: data_typecode = NPY_BOOL; break;
  case 2: data_typecode = NPY_UINT16; break;
  case 4: data_typecode = NPY_UINT32; break;
  case 8: data_typecode = NPY_UINT64; break;
  }
  return data_typecode;
}
/* -------------------------------------------------------------------------- */

template <>
inline int PythonFunctor::getPythonDataTypeCode<double>() const{
  return NPY_DOUBLE;
}

/* -------------------------------------------------------------------------- */

template <typename T>
PyObject *PythonFunctor::convertToNumpy(const T & akantu_object) const{
  AKANTU_DEBUG_TO_IMPLEMENT();
}
/* -------------------------------------------------------------------------- */

template <>
inline PyObject *PythonFunctor::convertToNumpy<double>(const double & akantu_object) const{
  return PyFloat_FromDouble(akantu_object);
}
/* -------------------------------------------------------------------------- */

template <>
inline PyObject *PythonFunctor::convertToNumpy<UInt>(const UInt & akantu_object) const{
  return PyInt_FromLong(akantu_object);
}
/* -------------------------------------------------------------------------- */

template <>
inline PyObject *PythonFunctor::convertToNumpy<bool>(const bool & akantu_object) const{
  return PyBool_FromLong(long(akantu_object));
}
/* -------------------------------------------------------------------------- */

template <typename T>
PyObject *PythonFunctor::convertToNumpy(const Vector<T> & array) const{
  int data_typecode = getPythonDataTypeCode< T >();
  npy_intp dims[1] = {array.size()};
  PyObject* obj = PyArray_SimpleNewFromData(1, dims, data_typecode, array.storage());
  PyArrayObject* res = (PyArrayObject*) obj;
  return (PyObject*)res;
}
/* -------------------------------------------------------------------------- */


inline PyObject * PythonFunctor::getPythonFunction(const std::string & functor_name) const{

  if (!PyInstance_Check(this->python_obj))
    AKANTU_EXCEPTION("Python object is not an instance");
  
  // if (!PyDict_Check(this->python_obj))
  //   AKANTU_EXCEPTION("Python object is not a dictionary");
  
  //PyObject * keys = PyDict_Keys(dict);
  //PyObject_Print(keys, stdout, Py_PRINT_RAW);    
  
  PyObject * pFunctor = PyObject_GetAttrString(this->python_obj, functor_name.c_str());
  if (!pFunctor) AKANTU_EXCEPTION("Python dictionary has no " << functor_name << " entry");

  return pFunctor;
}
/* -------------------------------------------------------------------------- */

inline void PythonFunctor::packArguments(std::vector<PyObject*> & pArgs) const{}

/* -------------------------------------------------------------------------- */


template<typename T, typename... Args>
inline void PythonFunctor::packArguments(std::vector<PyObject*> & pArgs,
					 T & p,
					 Args&...params) const{

  pArgs.push_back(this->convertToNumpy(p));
  if (sizeof...(params) != 0)
    this->packArguments(pArgs,params...);
}

/* -------------------------------------------------------------------------- */


template <typename return_type, typename... Params>
return_type PythonFunctor::callFunctor(const std::string & functor_name,
				       Params... parameters) const{


  _import_array();
  
  std::vector<PyObject*> arg_vector;
  this->packArguments(arg_vector,parameters...);

  PyObject * pArgs = PyTuple_New(arg_vector.size());
  for (UInt i = 0; i < arg_vector.size(); ++i) {
    PyTuple_SetItem(pArgs,i,arg_vector[i]);
  }
  
  PyObject * kwargs = PyDict_New();

  PyObject * pFunctor = getPythonFunction(functor_name);
  this->callFunctor(pFunctor,pArgs,kwargs);
  
  //Py_XDECREF(pArgs);
  //Py_XDECREF(kwargs);

}


__END_AKANTU__


#endif /* __AKANTU_PYTHON_FUNCTOR_INLINE_IMPL_CC__ */
