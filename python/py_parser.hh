#include <pybind11/pybind11.h>

#ifndef __AKANTU_PY_AKA_PARSER_HH__
#define __AKANTU_PY_AKA_PARSER_HH__

namespace akantu {

void register_parser(pybind11::module & mod);

}

#endif
