#ifndef __AKANTU_PYTHON_FUNCTOR_INLINE_IMPL_CC__
#define __AKANTU_PYTHON_FUNCTOR_INLINE_IMPL_CC__
/* -------------------------------------------------------------------------- */
#include <numpy/arrayobject.h>
#include "integration_point.hh"
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
PyObject *PythonFunctor::convertToPython(const T & akantu_object) const{
  AKANTU_DEBUG_TO_IMPLEMENT();
}
/* -------------------------------------------------------------------------- */

template <>
inline PyObject *PythonFunctor::convertToPython<double>(const double & akantu_object) const{
  return PyFloat_FromDouble(akantu_object);
}
/* -------------------------------------------------------------------------- */

template <>
inline PyObject *PythonFunctor::convertToPython<UInt>(const UInt & akantu_object) const{
  return PyInt_FromLong(akantu_object);
}
/* -------------------------------------------------------------------------- */

template <>
inline PyObject *PythonFunctor::convertToPython<bool>(const bool & akantu_object) const{
  return PyBool_FromLong(long(akantu_object));
}
/* -------------------------------------------------------------------------- */

template <typename T>
inline PyObject *PythonFunctor::convertToPython(const std::vector<T> & array) const{
  int data_typecode = getPythonDataTypeCode< T >();
  npy_intp dims[1] = {int(array.size())};
  PyObject* obj = PyArray_SimpleNewFromData(1, dims, data_typecode, const_cast<T*>(&array[0]));
  PyArrayObject* res = (PyArrayObject*) obj;
  return (PyObject*)res;
}
/* -------------------------------------------------------------------------- */

template <typename T>
PyObject *PythonFunctor::convertToPython(const Vector<T> & array) const{
  int data_typecode = getPythonDataTypeCode< T >();
  npy_intp dims[1] = {array.size()};
  PyObject* obj = PyArray_SimpleNewFromData(1, dims, data_typecode, array.storage());
  PyArrayObject* res = (PyArrayObject*) obj;
  return (PyObject*)res;
}
/* -------------------------------------------------------------------------- */

template <typename T>
PyObject *PythonFunctor::convertToPython(const Matrix<T> & mat) const{
  int data_typecode = getPythonDataTypeCode< T >();
  npy_intp dims[2] = {mat.size(0),mat.size(1)};
  PyObject* obj = PyArray_SimpleNewFromData(2, dims, data_typecode, mat.storage());
  PyArrayObject* res = (PyArrayObject*) obj;
  return (PyObject*)res;
}
/* -------------------------------------------------------------------------- */

template <>
inline PyObject *PythonFunctor::convertToPython<std::string>(const std::string & str) const{
  return PyString_FromString(str.c_str());  
}
/* -------------------------------------------------------------------------- */

template <>
inline PyObject *PythonFunctor::convertToPython<IntegrationPoint>(const IntegrationPoint & qp) const{

  PyObject * input = PyDict_New();
  PyObject * num_point   = this->convertToPython(qp.num_point);
  PyObject * global_num  = this->convertToPython(qp.global_num);
  PyObject * material_id = this->convertToPython(qp.material_id);
  PyObject * position    = this->convertToPython(qp.getPosition());
  PyDict_SetItemString(input,"num_point",num_point);
  PyDict_SetItemString(input,"global_num",global_num);
  PyDict_SetItemString(input,"material_id",material_id);
  PyDict_SetItemString(input,"position",position);
  return input;
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

inline void PythonFunctor::packArguments(std::vector<PyObject*> & pArgs)const {}

/* -------------------------------------------------------------------------- */


template<typename T, typename... Args>
inline void PythonFunctor::packArguments(std::vector<PyObject*> & pArgs,
					 T & p,
					 Args&...params)const{

  pArgs.push_back(this->convertToPython(p));
  if (sizeof...(params) != 0)
    this->packArguments(pArgs,params...);
}

/* -------------------------------------------------------------------------- */


template <typename return_type, typename... Params>
return_type PythonFunctor::callFunctor(const std::string & functor_name,
				       Params&... parameters) const{


  _import_array();
  
  std::vector<PyObject*> arg_vector;
  this->packArguments(arg_vector,parameters...);

  PyObject * pArgs = PyTuple_New(arg_vector.size());
  for (UInt i = 0; i < arg_vector.size(); ++i) {
    PyTuple_SetItem(pArgs,i,arg_vector[i]);
  }
  
  PyObject * kwargs = PyDict_New();

  PyObject * pFunctor = getPythonFunction(functor_name);
  PyObject * res = this->callFunctor(pFunctor,pArgs,kwargs);

  for (auto a: arg_vector) {
    // if (PyDict_Check(a)){
    //   PyObject* values = PyDict_Values(a);
    //   UInt sz = PyList_GET_SIZE(values);
    //   for (UInt i = 0; i < sz; ++i) {
    // 	Py_XDECREF(PyList_GetItem(values,i));
    //   }
    // }
    // Py_XDECREF(a);
  }
  Py_XDECREF(pArgs);
  Py_XDECREF(kwargs);

  return this->convertToAkantu<return_type>(res);

}

/* -------------------------------------------------------------------------- */

template <typename return_type>
inline return_type PythonFunctor::convertToAkantu(PyObject * python_obj) const{

  if (PyList_Check(python_obj)){
    return this->convertListToAkantu<typename return_type::value_type>(python_obj);
  }
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */

template <>
inline void PythonFunctor::convertToAkantu<void>(PyObject * python_obj) const{
  if (python_obj != Py_None) AKANTU_DEBUG_WARNING("functor return a value while none was expected: ignored");
}

/* -------------------------------------------------------------------------- */


template <>
inline std::string PythonFunctor::convertToAkantu<std::string>(PyObject * python_obj) const{
  if (!PyString_Check(python_obj)) AKANTU_EXCEPTION("cannot convert object to string");
  return PyString_AsString(python_obj);
}

/* -------------------------------------------------------------------------- */

template <>
inline Real PythonFunctor::convertToAkantu<Real>(PyObject * python_obj) const{
  if (!PyFloat_Check(python_obj)) AKANTU_EXCEPTION("cannot convert object to float");
  return PyFloat_AsDouble(python_obj);
}

/* -------------------------------------------------------------------------- */

template <typename T>
inline std::vector<T> PythonFunctor::convertListToAkantu(PyObject * python_obj) const{

  std::vector<T> res;
  UInt size = PyList_Size(python_obj);
  for (UInt i = 0; i < size; ++i) {
    PyObject * item = PyList_GET_ITEM(python_obj,i);
    res.push_back(this->convertToAkantu<T>(item));
  }
  return res;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__


#endif /* __AKANTU_PYTHON_FUNCTOR_INLINE_IMPL_CC__ */
