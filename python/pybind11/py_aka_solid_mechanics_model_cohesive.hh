#ifndef __AKANTU_PY_AKA_SOLID_MECHANICS_MODEL_COHESIVE_HH__
#define __AKANTU_PY_AKA_SOLID_MECHANICS_MODEL_COHESIVE_HH__

namespace pybind11 {
struct module;
}

namespace akantu {

void register_solid_mechanics_model_cohesive(pybind11::module & mod);

}

#endif
