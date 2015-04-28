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

#include <CGAL/Cartesian.h>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

EmbeddedInterfaceModel::EmbeddedInterfaceModel(Mesh & mesh,
                                               UInt spatial_dimension,
                                               const ID & id,
                                               const MemoryID & memory_id) :
  SolidMechanicsModel(mesh, spatial_dimension, id, memory_id),
  interface_mesh(NULL),
  interface_material_selector(NULL),
  interface_container(mesh)
{}

EmbeddedInterfaceModel::~EmbeddedInterfaceModel() {
  delete interface_material_selector;
}

void EmbeddedInterfaceModel::initModel() {
  interface_container.constructData();
  instanciateInterfaces();

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("reinforcement", id);
  this->mesh.addDumpMeshToDumper("reinforcement", *interface_mesh,
                                 1, _not_ghost, _ek_regular);
#endif

  SolidMechanicsModel::initModel();
}

void EmbeddedInterfaceModel::initMaterials() {
  delete interface_material_selector;
  interface_material_selector = new EmbeddedInterfaceMaterialSelector<std::string>("material", *this);
  this->setMaterialSelector(*interface_material_selector);

  Element element;
  //Material ** mat_val = &(materials.at(0));

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

void EmbeddedInterfaceModel::instanciateInterfaces() {
  const std::pair<Parser::const_section_iterator, Parser::const_section_iterator> &
    sub_sect = this->parser->getSubSections(_st_embedded_interface);

  Parser::const_section_iterator section_it = sub_sect.first;

  std::list<std::pair<K::Segment_3, std::string> > interface_list;

  for (; section_it != sub_sect.second ; ++section_it) {
    EmbeddedInterface interface(this->spatial_dimension);
    interface.parseSection(*section_it);

    interface_list.push_back(std::make_pair(interface.getPrimitive(), interface.getMaterialName()));
  }

  initInterface(interface_list);
}

void EmbeddedInterfaceModel::initInterface(const std::list<std::pair<K::Segment_3, std::string> > & interface_list) {
  AKANTU_DEBUG_IN();

  interface_mesh = &(interface_container.meshOfLinearInterfaces(interface_list));

  registerFEEngineObject<MyFEEngineType>("EmbeddedInterfaceFEEngine", *interface_mesh, 1);
  FEEngine & engine = getFEEngine("EmbeddedInterfaceFEEngine");
  engine.initShapeFunctions(_not_ghost);
  engine.initShapeFunctions(_ghost);

  AKANTU_DEBUG_OUT();
}

void EmbeddedInterfaceModel::assembleStiffnessMatrix() {
  SolidMechanicsModel::assembleStiffnessMatrix();
}

dumper::Field * EmbeddedInterfaceModel::createElementalField(const std::string & field_name,
    const std::string & group_name,
    bool padding_flag,
    const ElementKind & kind,
    const std::string & fe_engine_id) {

  if (field_name == "stress_embedded") {
    return SolidMechanicsModel::createElementalField(field_name,
                                                     group_name,
                                                     padding_flag,
                                                     kind,
                                                     "EmbeddedInterfaceFEEngine");
  } else {
    return SolidMechanicsModel::createElementalField(field_name,
                                                     group_name,
                                                     padding_flag,
                                                     kind,
                                                     fe_engine_id);
  }
}

ElementTypeMap<UInt> EmbeddedInterfaceModel::getInternalDataPerElem(const std::string & field_name,
                                                                    const ElementKind & kind,
                                                                    const std::string & fe_engine_id) {
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
        return materials[m]->getInternalDataPerElem(field_name,kind,fe_engine_id);
      else if (mat != NULL && field_name == "stress_embedded")
        return mat->getInternalDataPerElem(field_name, kind, fe_engine_id);
    }
  }

  return ElementTypeMap<UInt>();
}

__END_AKANTU__

