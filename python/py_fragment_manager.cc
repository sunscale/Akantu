/* -------------------------------------------------------------------------- */
#include <fragment_manager.hh>
#include <solid_mechanics_model_cohesive.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
#define def_deprecated(func_name, mesg)                                        \
  def(func_name, [](py::args, py::kwargs) { AKANTU_ERROR(mesg); })

#define def_function_nocopy(func_name)                                         \
  def(                                                                         \
      #func_name,                                                              \
      [](SolidMechanicsModel & self) -> decltype(auto) {                       \
        return self.func_name();                                               \
      },                                                                       \
      py::return_value_policy::reference)

#define def_function(func_name)                                                \
  def(#func_name, [](FragmentManager & self) -> decltype(auto) {               \
    return self.func_name();                                                   \
  })

void register_fragment_manager(py::module & mod) {
  py::class_<FragmentManager>(mod, "FragmentManager")
      .def(py::init<SolidMechanicsModelCohesive &, bool, const ID &>(),
           py::arg("model"), py::arg("dump_data") = true,
           py::arg("ID") = "fragment_manager")
      .def_function(buildFragments)
      .def_function(computeCenterOfMass)
      .def_function(computeVelocity)
      .def_function(computeInertiaMoments)
      .def_function(computeAllData)
      .def_function(computeNbElementsPerFragment)
      .def_function(getNbFragment)
      .def_function(getMass)
      .def_function(getVelocity)
      .def_function(getMomentsOfInertia)
      .def_function(getPrincipalDirections)
      .def_function(getNbElementsPerFragment);
}
} // namespace akantu
