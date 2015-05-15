/**
 * @file   embedded_interface_model.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Mon Mar 9 2015
 * @date last modification: Mon Mar 9 2015
 *
 * @brief  Model of Solid Mechanics with embedded interfaces
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#include "aka_common.hh"

#include "embedded_interface_model.hh"
#include "material_reinforcement.hh"

#ifdef AKANTU_USE_IOHELPER
#  include "dumper_paraview.hh"
#  include "dumpable_inline_impl.hh"
#endif

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

EmbeddedInterfaceModel::EmbeddedInterfaceModel(Mesh & mesh,
                                               Mesh & primitive_mesh,
                                               UInt spatial_dimension,
                                               const ID & id,
                                               const MemoryID & memory_id) :
  SolidMechanicsModel(mesh, spatial_dimension, id, memory_id),
  interface_mesh(NULL),
  primitive_mesh(primitive_mesh),
  interface_material_selector(NULL),
  intersector(mesh, primitive_mesh)
{
  // This pointer should be deleted by ~SolidMechanicsModel()
  MaterialSelector * mat_sel_pointer = new MeshDataMaterialSelector<std::string>("physical_names", *this);
  this->setMaterialSelector(*mat_sel_pointer);
}

EmbeddedInterfaceModel::~EmbeddedInterfaceModel() {
  delete interface_material_selector;
}

void EmbeddedInterfaceModel::initFull(const ModelOptions & options) {
  const EmbeddedInterfaceModelOptions & eim_options =
    dynamic_cast<const EmbeddedInterfaceModelOptions &>(options);

  // We don't want to initiate materials before shape functions are initialized
  SolidMechanicsModelOptions dummy_options(eim_options.analysis_method, true);

  interface_mesh = &(intersector.getInterfaceMesh());
  registerFEEngineObject<MyFEEngineType>("EmbeddedInterfaceFEEngine", *interface_mesh, 1);

  SolidMechanicsModel::initFull(dummy_options);

  if (!eim_options.no_init_intersections)
    intersector.constructData();

  FEEngine & engine = getFEEngine("EmbeddedInterfaceFEEngine");
  engine.initShapeFunctions(_not_ghost);
  engine.initShapeFunctions(_ghost);

  this->initMaterials();

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("reinforcement", id);
  this->mesh.addDumpMeshToDumper("reinforcement", *interface_mesh,
                                 1, _not_ghost, _ek_regular);
#endif
}

void EmbeddedInterfaceModel::initMaterials() {
  Element element;

  delete interface_material_selector;
  interface_material_selector = new InterfaceMeshDataMaterialSelector<std::string>("physical_names", *this);

  for (ghost_type_t::iterator gt = ghost_type_t::begin(); gt != ghost_type_t::end(); ++gt) {
    element.ghost_type = *gt;

    Mesh::type_iterator it = interface_mesh->firstType(1, *gt);
    Mesh::type_iterator end = interface_mesh->lastType(1, *gt);

    for (; it != end ; ++it) {
      UInt nb_element = interface_mesh->getNbElement(*it, *gt);

      element.type = *it;

      Array<UInt> & mat_indexes = material_index.alloc(nb_element, 1, *it, *gt);

      for (UInt el = 0 ; el < nb_element ; el++) {
        element.element = el;
        UInt mat_index = (*interface_material_selector)(element);

        AKANTU_DEBUG_ASSERT(mat_index < materials.size(),
            "The material selector returned an index that does not exist");
        mat_indexes(element.element) = mat_index;
        materials.at(mat_index)->addElement(*it, el, *gt);
      }
    }
  }

  SolidMechanicsModel::initMaterials();
}

ElementTypeMap<UInt> EmbeddedInterfaceModel::getInternalDataPerElem(const std::string & field_name,
                                                                    const ElementKind & kind) {
  if (!(this->isInternal(field_name,kind))) AKANTU_EXCEPTION("unknown internal " << field_name);

  for (UInt m = 0; m < materials.size() ; ++m) {
    if (materials[m]->isInternal(field_name, kind)) {
      Material * mat = NULL;

      switch(this->spatial_dimension) {
        case 1:
          mat = dynamic_cast<MaterialReinforcement<1> *>(materials[m]);
          break;

        case 2:
          mat = dynamic_cast<MaterialReinforcement<2> *>(materials[m]);
          break;

        case 3:
          mat = dynamic_cast<MaterialReinforcement<3> *>(materials[m]);
          break;
      }

      if (mat == NULL && field_name != "stress_embedded")
        return materials[m]->getInternalDataPerElem(field_name,kind);
      else if (mat != NULL && field_name == "stress_embedded")
        return mat->getInternalDataPerElem(field_name, kind, "EmbeddedInterfaceFEEngine");
    }
  }

  return ElementTypeMap<UInt>();
}

void EmbeddedInterfaceModel::addDumpGroupFieldToDumper(const std::string & dumper_name,
                                                       const std::string & field_id,
                                                       const std::string & group_name,
                                                       const ElementKind & element_kind,
                                                       bool padding_flag) {
#ifdef AKANTU_USE_IOHELPER
  dumper::Field * field = NULL;

  if (field_id == "stress_embedded")
    field = this->createElementalField(field_id, group_name, padding_flag, 1, element_kind);
  else
    SolidMechanicsModel::addDumpGroupFieldToDumper(dumper_name, field_id, group_name, element_kind, padding_flag);

  if (field) {
    DumperIOHelper & dumper = mesh.getGroupDumper(dumper_name,group_name);
    Model::addDumpGroupFieldToDumper(field_id,field,dumper);
  }

#endif
}

__END_AKANTU__

