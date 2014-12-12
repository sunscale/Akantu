%module akantu

%include "stl.i"

#define __attribute__(x)

%ignore akantu::operator <<;

%include "aka_common.i"
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

