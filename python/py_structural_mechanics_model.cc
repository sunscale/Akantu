/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include "structural_mechanics_model.hh"
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
	def(#func_name,							       \
	    [](StructuralMechanicsModel & self) -> decltype(auto)       	       \
				{return self.func_name();},		       \
	    py::return_value_policy::reference)

#define def_function_(func_name)                                               \
	def(#func_name,							       \
	    [](StructuralMechanicsModel & self) -> decltype(auto)       	       \
				{return self.func_name();})

#define def_plainmember(M) def_readwrite(#M, &StructuralMaterial:: M)
/* -------------------------------------------------------------------------- */


void
register_structural_mechanics_model(
	pybind11::module & mod)
{
	/* First we have to register the material class
	 *  The wrapper aims to mimic the behaviour of the real material.
	 */
	py::class_<StructuralMaterial>(mod, "StructuralMaterial")
		.def_plainmember(E)
		.def_plainmember(A)
		.def_plainmember(I)
		.def_plainmember(Iz)
		.def_plainmember(Iy)
		.def_plainmember(GJ)
		.def_plainmember(rho)
		.def_plainmember(t)
		.def_plainmember(nu);


	/* Now we create the structural model wrapper
	 *  Note that this is basically a port from the solid mechanic part.
	 */
	py::class_<StructuralMechanicsModel, Model>(mod, "StructuralMechanicsModel")
		.def(py::init<Mesh &, UInt, const ID &, const MemoryID &, const ModelType>(),
			py::arg("mesh"),
			py::arg("spatial_dimension") = _all_dimensions,
			py::arg("id") = "structural_mechanics_model",
			py::arg("memory_id") = 0
		)
		.def("initFull",
			[](StructuralMechanicsModel& self, const AnalysisMethod& analysis_method) -> void
			   { self.initFull(_analysis_method = analysis_method);	return; },
			py::arg("_analysis_method")
		)
		.def("setTimeStep",
			[](StructuralMechanicsModel& self, const Real& time_step, const ID& solver_id) -> void
			   { AKANTU_ERROR("This function was commented out in the source code."); return;
			     (void)self; (void)time_step; (void)solver_id; },
			py::arg("time_step"),
			py::arg("solver_id") = ""
		)
		.def_function_nocopy(getExternalForce)
		.def_function_nocopy(getDisplacement)
		.def_function_nocopy(getInternalForce)
		.def_function_nocopy(getVelocity)
		.def_function_nocopy(getAcceleration)
		.def_function_nocopy(getInternalForce)
		.def_function_nocopy(getBlockedDOFs)
		.def_function_nocopy(getMesh)
#if 0
		.def("dump", py::overload_cast<>(&StructuralMechanicsModel::dump))
		.def("dump",
				py::overload_cast<const std::string &>(&StructuralMechanicsModel::dump))
		.def("dump", py::overload_cast<const std::string &, UInt>(
					&StructuralMechanicsModel::dump))
		.def("dump", py::overload_cast<const std::string &, Real, UInt>(
					&StructuralMechanicsModel::dump))
		.def("getMaterial",
				py::overload_cast<UInt>(&StructuralMechanicsModel::getMaterial),
				py::return_value_policy::reference)
		.def("getMaterial",
				py::overload_cast<const std::string &>(
					&StructuralMechanicsModel::getMaterial),
				py::return_value_policy::reference)
		.def("getMaterialIndex", &StructuralMechanicsModel::getMaterialIndex)
		.def("setMaterialSelector", &StructuralMechanicsModel::setMaterialSelector)
		.def("getMaterialSelector", &StructuralMechanicsModel::getMaterialSelector)
#endif
		;

}; //End: register structural mechanical model

} // namespace akantu

