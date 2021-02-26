/**
 * @file   embedded_interface_model.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Mar 13 2015
 * @date last modification: Wed Feb 14 2018
 *
 * @brief  Model of Solid Mechanics with embedded interfaces
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#include "embedded_interface_model.hh"
#include "integrator_gauss.hh"
#include "material_elastic.hh"
#include "material_reinforcement.hh"
#include "mesh_iterators.hh"
#include "shape_lagrange.hh"

#ifdef AKANTU_USE_IOHELPER
#include "dumpable_inline_impl.hh"
#include "dumper_iohelper_paraview.hh"
#endif

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
EmbeddedInterfaceModel::EmbeddedInterfaceModel(Mesh & mesh,
                                               Mesh & primitive_mesh,
                                               UInt spatial_dimension,
                                               const ID & id,
                                               const MemoryID & memory_id)
    : SolidMechanicsModel(mesh, spatial_dimension, id, memory_id),
      intersector(mesh, primitive_mesh), interface_mesh(nullptr),
      primitive_mesh(primitive_mesh), interface_material_selector(nullptr) {
  this->model_type = ModelType::_embedded_model;

  // This pointer should be deleted by ~SolidMechanicsModel()
  auto mat_sel_pointer =
      std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names",
                                                              *this);

  this->setMaterialSelector(mat_sel_pointer);

  interface_mesh = &(intersector.getInterfaceMesh());

  // Create 1D FEEngine on the interface mesh
  registerFEEngineObject<MyFEEngineType>("EmbeddedInterfaceFEEngine",
                                         *interface_mesh, 1);

  // Registering allocator for material reinforcement
  MaterialFactory::getInstance().registerAllocator(
      "reinforcement",
      [&](UInt dim, const ID & constitutive, SolidMechanicsModel & /*unused*/,
          const ID & id) -> std::unique_ptr<Material> {
        if (constitutive == "elastic") {
          using mat = MaterialElastic<1>;
          switch (dim) {
          case 2:
            return std::make_unique<MaterialReinforcement<mat, 2>>(*this, id);
          case 3:
            return std::make_unique<MaterialReinforcement<mat, 3>>(*this, id);
          default:
            AKANTU_EXCEPTION("Dimension 1 is invalid for reinforcements");
          }
        } else {
          AKANTU_EXCEPTION("Reinforcement type" << constitutive
                                                << " is not recognized");
        }
      });
}

/* -------------------------------------------------------------------------- */
EmbeddedInterfaceModel::~EmbeddedInterfaceModel() {
  delete interface_material_selector;
}

/* -------------------------------------------------------------------------- */
void EmbeddedInterfaceModel::initFullImpl(const ModelOptions & options) {
  const auto & eim_options =
      aka::as_type<EmbeddedInterfaceModelOptions>(options);

  // Do no initialize interface_mesh if told so
  if (eim_options.has_intersections) {
    intersector.constructData();
  }

  SolidMechanicsModel::initFullImpl(options);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("reinforcement", id);
  this->mesh.addDumpMeshToDumper("reinforcement", *interface_mesh, 1,
                                 _not_ghost, _ek_regular);
#endif
}

void EmbeddedInterfaceModel::initModel() {
  // Initialize interface FEEngine
  SolidMechanicsModel::initModel();
  FEEngine & engine = getFEEngine("EmbeddedInterfaceFEEngine");
  engine.initShapeFunctions(_not_ghost);
  engine.initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */
void EmbeddedInterfaceModel::assignMaterialToElements(
    const ElementTypeMapArray<UInt> * filter) {
  delete interface_material_selector;
  interface_material_selector =
      new InterfaceMeshDataMaterialSelector<std::string>("physical_names",
                                                         *this);

  for_each_element(getInterfaceMesh(),
                   [&](auto && element) {
                     auto mat_index = (*interface_material_selector)(element);
                     // material_index(element) = mat_index;
                     materials[mat_index]->addElement(element);
                     // this->material_local_numbering(element) = index;
                   },
                   _element_filter = filter, _spatial_dimension = 1);

  SolidMechanicsModel::assignMaterialToElements(filter);
}

/* -------------------------------------------------------------------------- */
void EmbeddedInterfaceModel::addDumpGroupFieldToDumper(
    const std::string & dumper_name, const std::string & field_id,
    const std::string & group_name, ElementKind element_kind,
    bool padding_flag) {
#ifdef AKANTU_USE_IOHELPER
  std::shared_ptr<dumpers::Field> field;

  // If dumper is reinforcement, create a 1D elemental field
  if (dumper_name == "reinforcement") {
    field = this->createElementalField(field_id, group_name, padding_flag, 1,
                                       element_kind);
  } else {
    try {
      SolidMechanicsModel::addDumpGroupFieldToDumper(
          dumper_name, field_id, group_name, element_kind, padding_flag);
    } catch (...) {
    }
  }
  if (field) {
    DumperIOHelper & dumper = mesh.getGroupDumper(dumper_name, group_name);
    Model::addDumpGroupFieldToDumper(field_id, field, dumper);
  }

#endif
}

} // namespace akantu
