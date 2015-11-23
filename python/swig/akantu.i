%module akantu

%exception {
  try {
    $action
  } catch (akantu::debug::Exception e) {
    PyErr_SetString(PyExc_IndexError,e.what());
    return NULL;
  }
 }

%include "stl.i"

#define __attribute__(x)

%ignore akantu::operator <<;

%include "aka_common.i"
%include "aka_csr.i"
%include "aka_array.i"

%define print_self(MY_CLASS)
  %extend akantu::MY_CLASS {
    std::string __str__() {
      std::stringstream sstr;
      sstr << *($self);
      return sstr.str();
    }
 }
%enddef

%include "mesh.i"
%include "mesh_utils.i"
%include "model.i"
%include "solid_mechanics_model.i"
#if defined(AKANTU_COHESIVE_ELEMENT)
%include "solid_mechanics_model_cohesive.i"
#endif

#if defined(AKANTU_HEAT_TRANSFER)
%include "heat_transfer_model.i"
#endif


#if defined(AKANTU_STRUCTURAL_MECHANICS)
%include "load_functions.i"
%include "structural_mechanics_model.i"
#endif
