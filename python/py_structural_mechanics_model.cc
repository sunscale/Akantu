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
	    [](StructuralMechanicsModel & self) -> decltype(auto)              \
				{return self.func_name();},		       \
	    py::return_value_policy::reference)

#define def_function_(func_name)                                               \
	def(#func_name,							       \
	    [](StructuralMechanicsModel & self) -> decltype(auto)              \
				{return self.func_name();})

#define def_plainmember(M) 						       \
	def_readwrite(#M, &StructuralMaterial:: M)
/* -------------------------------------------------------------------------- */


void
register_structural_mechanics_model(
	pybind11::module & mod)
{
	/* First we have to register the material class
	 *  The wrapper aims to mimic the behaviour of the real material.
	 */
	py::class_<StructuralMaterial>(mod, "StructuralMaterial")
		.def(py::init<>())
		.def_plainmember(E)
		.def_plainmember(A)
		.def_plainmember(I)
		.def_plainmember(Iz)
		.def_plainmember(Iy)
		.def_plainmember(GJ)
		.def_plainmember(rho)
		.def_plainmember(t)
		.def_plainmember(nu)
	;


	/* Now we create the structural model wrapper
	 *  Note that this is basically a port from the solid mechanic part.
	 */
	py::class_<StructuralMechanicsModel, Model>(mod, "StructuralMechanicsModel")
		.def(py::init<Mesh &, UInt, const ID &, const MemoryID &>(),
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
		.def("initFull",
			[](StructuralMechanicsModel& self) -> void
			   { self.initFull();	return; }
		)
#		if 0
		.def("setTimeStep",
			[](StructuralMechanicsModel& self, const Real& time_step, const ID& solver_id) -> void
			   { AKANTU_ERROR("This function was commented out in the source code."); return;
			     (void)self; (void)time_step; (void)solver_id; },
			py::arg("time_step"),
			py::arg("solver_id") = ""
		)
#		endif
		.def_function_nocopy(getExternalForce)
		.def_function_nocopy(getDisplacement)
		.def_function_nocopy(getInternalForce)
		.def_function_nocopy(getVelocity)
		.def_function_nocopy(getAcceleration)
		.def_function_nocopy(getInternalForce)
		.def_function_nocopy(getBlockedDOFs)
		.def_function_nocopy(getMesh)

		/*
		 * These functions are basically untested.
		 */
		.def("getElementMaterialMap",
			[](StructuralMechanicsModel& self, const ElementType& type, GhostType ghost_type)
			   { return self.getElementMaterial(type, ghost_type); },
			"This function returns the map that maps elements to materials.",
			py::arg("type"),
			py::arg("ghost_type") = _not_ghost,
          		py::return_value_policy::reference
		)
		.def("getMaterialOf",
			[](StructuralMechanicsModel& self, Element element)
			   { return self.getMaterial(element); },
			"This function returns the `StructuralMaterial` instance that is associated with element `element`."
			" It is important that the returned object can be modified, but this will not affect the material stored inside the model."
			" If you want to change the material, use `addMaterial()` to add a new one and then manipulate the mapping by operating on `getElementMaterialMap()`.",
			py::arg("element"),
			py::return_value_policy::copy	//By using the copy operation, we completly decouple the C++ and Python part.
		)
		.def("addMaterial",
			[](StructuralMechanicsModel& self, StructuralMaterial& mat) -> UInt
			   { return self.addMaterial(mat); },
			"This function adds the `StructuralMaterial` `mat` to `self`."
			" The function returns the ID of the new material.",
			py::arg("mat")
		)
		.def("getMaterialByID",
			[](StructuralMechanicsModel& self, UInt i)
			   { return self.getMaterialByID(i); },
			"This function returns the `i`th material of `self`",
			py::arg("i"),
			py::return_value_policy::copy 	//Everything will be coupled so a complet decoupling
		)
		.def("getNbMaterials",
			[](StructuralMechanicsModel& self)
			   { return self.getNbMaterials(); },
			"Returns the number of different materials inside `self`."
			" The highest ID is one less than the returned number."
		)
		;

	return;
} //End: register structural mechanical model

} // namespace akantu

