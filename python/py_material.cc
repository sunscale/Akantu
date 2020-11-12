/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <solid_mechanics_model.hh>
#if defined(AKANTU_COHESIVE_ELEMENT)
#include <solid_mechanics_model_cohesive.hh>
#endif
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

#if not defined(PYBIND11_OVERRIDE)
#define PYBIND11_OVERRIDE PYBIND11_OVERLOAD
#endif

namespace akantu {

template <typename _Material> class PyMaterial : public _Material {

public:
  /* Inherit the constructors */
  using _Material::_Material;

  ~PyMaterial() override = default;

  void initMaterial() override {
    PYBIND11_OVERRIDE(void, _Material, initMaterial); // NOLINT
  };
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override {
    PYBIND11_OVERRIDE_PURE(void, _Material, computeStress, el_type, ghost_type);
  }
  void computeTangentModuli(ElementType el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override {
    PYBIND11_OVERRIDE(void, _Material, computeTangentModuli, el_type,
                      tangent_matrix, ghost_type);
  }

  void computePotentialEnergy(ElementType el_type) override {
    PYBIND11_OVERRIDE(void, _Material, computePotentialEnergy, el_type);
  }

  Real getPushWaveSpeed(const Element & element) const override {
    PYBIND11_OVERRIDE(Real, _Material, getPushWaveSpeed, element);
  }

  Real getShearWaveSpeed(const Element & element) const override {
    PYBIND11_OVERRIDE(Real, _Material, getShearWaveSpeed, element);
  }

  void registerInternal(const std::string & name, UInt nb_component) {
    this->internals[name] = std::make_shared<InternalField<Real>>(name, *this);
    AKANTU_DEBUG_INFO("alloc internal " << name << " "
                                        << &this->internals[name]);

    this->internals[name]->initialize(nb_component);
  }

  auto & getInternals() { return this->internals; }

protected:
  std::map<std::string, std::shared_ptr<InternalField<Real>>> internals;
};

/* -------------------------------------------------------------------------- */
template <typename T>
void register_element_type_map_array(py::module & mod,
                                     const std::string & name) {
  py::class_<ElementTypeMapArray<T>, std::shared_ptr<ElementTypeMapArray<T>>>(
      mod, ("ElementTypeMapArray" + name).c_str())
      .def(
          "__call__",
          [](ElementTypeMapArray<T> & self, ElementType type,
             GhostType ghost_type) -> decltype(auto) {
            return self(type, ghost_type);
          },
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "elementTypes",
          [](ElementTypeMapArray<T> & self, UInt _dim, GhostType _ghost_type,
             ElementKind _kind) -> decltype(auto) {
            auto types = self.elementTypes(_dim, _ghost_type, _kind);
            std::vector<ElementType> _types;
            for (auto && t : types) {
              _types.push_back(t);
            }
            return _types;
          },
          py::arg("dim") = _all_dimensions, py::arg("ghost_type") = _not_ghost,
          py::arg("kind") = _ek_regular);

  py::class_<InternalField<T>, ElementTypeMapArray<T>,
             std::shared_ptr<InternalField<T>>>(
      mod, ("InternalField" + name).c_str());
}

/* -------------------------------------------------------------------------- */
template <typename _Material>
void define_material(py::module & mod, const std::string & name) {
  py::class_<_Material, PyMaterial<_Material>, Parsable>(
      mod, name.c_str(), py::multiple_inheritance())
      .def(py::init<SolidMechanicsModel &, const ID &>())
      .def(
          "getGradU",
          [](Material & self, ElementType el_type,
             GhostType ghost_type = _not_ghost) -> decltype(auto) {
            return self.getGradU(el_type, ghost_type);
          },
          py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "getStress",
          [](Material & self, ElementType el_type,
             GhostType ghost_type = _not_ghost) -> decltype(auto) {
            return self.getStress(el_type, ghost_type);
          },
          py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "getPotentialEnergy",
          [](Material & self, ElementType el_type) -> decltype(auto) {
            return self.getPotentialEnergy(el_type);
          },
          py::return_value_policy::reference)
      .def("initMaterial", &Material::initMaterial)
      .def("getModel", &Material::getModel)
      .def("registerInternal",
           [](Material & self, const std::string & name, UInt nb_component) {
             return dynamic_cast<PyMaterial<Material> &>(self).registerInternal(
                 name, nb_component);
           })
      .def(
          "getInternalFieldReal",
          [](Material & self, const ID & id, const ElementType & type,
             const GhostType & ghost_type) -> Array<Real> & {
            return self.getArray<Real>(id, type, ghost_type);
          },
          py::arg("id"), py::arg("type"), py::arg("ghost_type") = _not_ghost)
      .def(
          "getInternalFieldUInt",
          [](Material & self, const ID & id, const ElementType & type,
             const GhostType & ghost_type) -> Array<UInt> & {
            return self.getArray<UInt>(id, type, ghost_type);
          },
          py::arg("id"), py::arg("type"), py::arg("ghost_type") = _not_ghost)
      .def(
          "getElementFilter",
          [](Material & self, const ElementType & type,
             const GhostType & ghost_type) -> const Array<UInt> & {
            return self.getElementFilter()(type, ghost_type);
          },
          py::arg("type"), py::arg("ghost_type") = _not_ghost);
}

/* -------------------------------------------------------------------------- */
void register_material(py::module & mod) {
  py::class_<MaterialFactory>(mod, "MaterialFactory")
      .def_static(
          "getInstance",
          []() -> MaterialFactory & { return Material::getFactory(); },
          py::return_value_policy::reference)
      .def("registerAllocator",
           [](MaterialFactory & self, const std::string id, py::function func) {
             self.registerAllocator(
                 id,
                 [func, id](UInt dim, const ID & /*unused*/,
                            SolidMechanicsModel & model,
                            const ID & option) -> std::unique_ptr<Material> {
                   py::object obj = func(dim, id, model, option);
                   auto & ptr = py::cast<Material &>(obj);

                   obj.release();
                   return std::unique_ptr<Material>(&ptr);
                 });
           });

  register_element_type_map_array<Real>(mod, "Real");
  register_element_type_map_array<UInt>(mod, "UInt");

  define_material<Material>(mod, "Material");

  py::class_<MaterialSelector, std::shared_ptr<MaterialSelector>>(
      mod, "MaterialSelector")
      .def("setFallback",
           [](MaterialSelector & self, UInt f) { self.setFallback(f); })
      .def("setFallback",
           [](MaterialSelector & self,
              const std::shared_ptr<MaterialSelector> & fallback_selector) {
             self.setFallback(fallback_selector);
           })
      .def("setFallback",
           [](MaterialSelector & self, MaterialSelector & fallback_selector) {
             self.setFallback(fallback_selector);
           });

  py::class_<MeshDataMaterialSelector<std::string>, MaterialSelector,
             std::shared_ptr<MeshDataMaterialSelector<std::string>>>(
      mod, "MeshDataMaterialSelectorString")
      .def(py::init<const std::string &, const SolidMechanicsModel &, UInt>(),
           py::arg("name"), py::arg("model"), py::arg("first_index") = 1);

#if defined(AKANTU_COHESIVE_ELEMENT)
  py::class_<DefaultMaterialCohesiveSelector, MaterialSelector,
             std::shared_ptr<DefaultMaterialCohesiveSelector>>(
      mod, "DefaultMaterialCohesiveSelector")
      .def(py::init<const SolidMechanicsModelCohesive &>());

  py::class_<MeshDataMaterialCohesiveSelector, MaterialSelector,
             std::shared_ptr<MeshDataMaterialCohesiveSelector>>(
      mod, "MeshDataMaterialCohesiveSelector")
      .def(py::init<const SolidMechanicsModelCohesive &>());

  py::class_<MaterialCohesiveRulesSelector, MaterialSelector,
             std::shared_ptr<MaterialCohesiveRulesSelector>>(
      mod, "MaterialCohesiveRulesSelector")
      .def(py::init<const SolidMechanicsModelCohesive &,
                    const MaterialCohesiveRules &, const ID &>(),
           py::arg("model"), py::arg("rules"),
           py::arg("mesh_data_id") = "physical_names");
#endif
}

} // namespace akantu
