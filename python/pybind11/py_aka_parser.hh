#ifndef __AKANTU_PY_AKA_PARSER_HH__
#define __AKANTU_PY_AKA_PARSER_HH__

namespace pybind11 {
struct module;
}

namespace akantu {

void register_parser(pybind11::module & mod);

}

#endif
