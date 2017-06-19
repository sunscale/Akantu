/**
 * @file   model.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Oct 03 2011
 * @date last modification: Thu Nov 19 2015
 *
 * @brief  implementation of model common parts
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "model.hh"
#include "element_group.hh"
#include "synchronizer_registry.hh"
#include "data_accessor.hh"
#include "static_communicator.hh"
#include "element_synchronizer.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Model::Model(Mesh& mesh, UInt dim, const ID & id,
             const MemoryID & memory_id) :
  Memory(id, memory_id), ModelSolver(mesh, id, memory_id),
  mesh(mesh),
  spatial_dimension(dim == _all_dimensions ? mesh.getSpatialDimension() : dim),
  is_pbc_slave_node(0,1,"is_pbc_slave_node") ,
  parser(&getStaticParser()) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Model::~Model() {
  AKANTU_DEBUG_IN();

  FEEngineMap::iterator it;
  for (it = fems.begin(); it != fems.end(); ++it) {
    if (it->second)
      delete it->second;
  }

  for (it = fems_boundary.begin(); it != fems_boundary.end(); ++it) {
    if (it->second)
      delete it->second;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Model::setParser(Parser & parser) { this->parser = &parser; }

/* -------------------------------------------------------------------------- */
void Model::initFull(__attribute__((unused)) const ModelOptions & options) {
  AKANTU_DEBUG_IN();

  initModel();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Model::initPBC() {
  std::map<UInt, UInt>::iterator it = pbc_pair.begin();
  std::map<UInt, UInt>::iterator end = pbc_pair.end();

  is_pbc_slave_node.resize(mesh.getNbNodes());
#ifndef AKANTU_NDEBUG
  Array<Real>::const_vector_iterator coord_it =
      mesh.getNodes().begin(this->spatial_dimension);
#endif

  while (it != end) {
    UInt i1 = (*it).first;

    is_pbc_slave_node(i1) = true;

#ifndef AKANTU_NDEBUG
    UInt i2 = (*it).second;
    UInt slave = mesh.isDistributed() ? mesh.getGlobalNodesIds()(i1) : i1;
    UInt master = mesh.isDistributed() ? mesh.getGlobalNodesIds()(i2) : i2;

    AKANTU_DEBUG_INFO("pairing " << slave << " (" << Vector<Real>(coord_it[i1])
                                 << ") with " << master << " ("
                                 << Vector<Real>(coord_it[i2]) << ")");
#endif
    ++it;
  }
}

/* -------------------------------------------------------------------------- */
void Model::dumpGroup(const std::string & group_name) {
  ElementGroup & group = mesh.getElementGroup(group_name);
  group.dump();
}

/* -------------------------------------------------------------------------- */
void Model::dumpGroup(const std::string & group_name,
                      const std::string & dumper_name) {
  ElementGroup & group = mesh.getElementGroup(group_name);
  group.dump(dumper_name);
}

/* -------------------------------------------------------------------------- */
void Model::dumpGroup() {
  GroupManager::element_group_iterator bit = mesh.element_group_begin();
  GroupManager::element_group_iterator bend = mesh.element_group_end();
  for (; bit != bend; ++bit) {
    bit->second->dump();
  }
}

/* -------------------------------------------------------------------------- */
void Model::setGroupDirectory(const std::string & directory) {
  GroupManager::element_group_iterator bit = mesh.element_group_begin();
  GroupManager::element_group_iterator bend = mesh.element_group_end();
  for (; bit != bend; ++bit) {
    bit->second->setDirectory(directory);
  }
}

/* -------------------------------------------------------------------------- */
void Model::setGroupDirectory(const std::string & directory,
                              const std::string & group_name) {
  ElementGroup & group = mesh.getElementGroup(group_name);
  group.setDirectory(directory);
}

/* -------------------------------------------------------------------------- */
void Model::setGroupBaseName(const std::string & basename,
                             const std::string & group_name) {
  ElementGroup & group = mesh.getElementGroup(group_name);
  group.setBaseName(basename);
}

/* -------------------------------------------------------------------------- */
DumperIOHelper & Model::getGroupDumper(const std::string & group_name) {
  ElementGroup & group = mesh.getElementGroup(group_name);
  return group.getDumper();
}

/* -------------------------------------------------------------------------- */
// DUMPER stuff
/* -------------------------------------------------------------------------- */
void Model::addDumpGroupFieldToDumper(const std::string & field_id,
                                      dumper::Field * field,
                                      DumperIOHelper & dumper) {
#ifdef AKANTU_USE_IOHELPER
  dumper.registerField(field_id, field);
#endif
}

/* -------------------------------------------------------------------------- */
void Model::addDumpField(const std::string & field_id) {

  this->addDumpFieldToDumper(mesh.getDefaultDumperName(), field_id);
}
/* -------------------------------------------------------------------------- */

void Model::addDumpFieldVector(const std::string & field_id) {

  this->addDumpFieldVectorToDumper(mesh.getDefaultDumperName(), field_id);
}

/* -------------------------------------------------------------------------- */
void Model::addDumpFieldTensor(const std::string & field_id) {

  this->addDumpFieldTensorToDumper(mesh.getDefaultDumperName(), field_id);
}

/* -------------------------------------------------------------------------- */

void Model::setBaseName(const std::string & field_id) {

  mesh.setBaseName(field_id);
}
/* -------------------------------------------------------------------------- */

void Model::setBaseNameToDumper(const std::string & dumper_name,
                                const std::string & basename) {
  mesh.setBaseNameToDumper(dumper_name, basename);
}
/* -------------------------------------------------------------------------- */

void Model::addDumpFieldToDumper(const std::string & dumper_name,
                                 const std::string & field_id) {

  this->addDumpGroupFieldToDumper(dumper_name, field_id, "all", _ek_regular,
                                  false);
}

/* -------------------------------------------------------------------------- */
void Model::addDumpGroupField(const std::string & field_id,
                              const std::string & group_name) {

  ElementGroup & group = mesh.getElementGroup(group_name);
  this->addDumpGroupFieldToDumper(group.getDefaultDumperName(), field_id,
                                  group_name, _ek_regular, false);
}

/* -------------------------------------------------------------------------- */
void Model::removeDumpGroupField(const std::string & field_id,
                                 const std::string & group_name) {
  ElementGroup & group = mesh.getElementGroup(group_name);
  this->removeDumpGroupFieldFromDumper(group.getDefaultDumperName(), field_id,
                                       group_name);
}

/* -------------------------------------------------------------------------- */
void Model::removeDumpGroupFieldFromDumper(const std::string & dumper_name,
                                           const std::string & field_id,
                                           const std::string & group_name) {
  ElementGroup & group = mesh.getElementGroup(group_name);
  group.removeDumpFieldFromDumper(dumper_name, field_id);
}

/* -------------------------------------------------------------------------- */
void Model::addDumpFieldVectorToDumper(const std::string & dumper_name,
                                       const std::string & field_id) {
  this->addDumpGroupFieldToDumper(dumper_name, field_id, "all", _ek_regular, 3);
}

/* -------------------------------------------------------------------------- */
void Model::addDumpGroupFieldVector(const std::string & field_id,
                                    const std::string & group_name) {
  ElementGroup & group = mesh.getElementGroup(group_name);
  this->addDumpGroupFieldVectorToDumper(group.getDefaultDumperName(), field_id,
                                        group_name);
}

/* -------------------------------------------------------------------------- */
void Model::addDumpGroupFieldVectorToDumper(const std::string & dumper_name,
                                            const std::string & field_id,
                                            const std::string & group_name) {
  this->addDumpGroupFieldToDumper(dumper_name, field_id, group_name,
                                  _ek_regular, true);
}
/* -------------------------------------------------------------------------- */

void Model::addDumpFieldTensorToDumper(const std::string & dumper_name,
                                       const std::string & field_id) {
  this->addDumpGroupFieldToDumper(dumper_name, field_id, "all", _ek_regular,
                                  true);
}

/* -------------------------------------------------------------------------- */
void Model::addDumpGroupFieldToDumper(const std::string & dumper_name,
                                      const std::string & field_id,
                                      const std::string & group_name,
                                      const ElementKind & element_kind,
                                      bool padding_flag) {
  this->addDumpGroupFieldToDumper(dumper_name, field_id, group_name,
                                  this->spatial_dimension, element_kind,
                                  padding_flag);
}

/* -------------------------------------------------------------------------- */
void Model::addDumpGroupFieldToDumper(const std::string & dumper_name,
                                      const std::string & field_id,
                                      const std::string & group_name,
                                      UInt spatial_dimension,
                                      const ElementKind & element_kind,
                                      bool padding_flag) {

#ifdef AKANTU_USE_IOHELPER
  dumper::Field * field = NULL;

  if (!field)
    field = this->createNodalFieldReal(field_id, group_name, padding_flag);
  if (!field)
    field = this->createNodalFieldUInt(field_id, group_name, padding_flag);
  if (!field)
    field = this->createNodalFieldBool(field_id, group_name, padding_flag);
  if (!field)
    field = this->createElementalField(field_id, group_name, padding_flag,
                                       spatial_dimension, element_kind);
  if (!field)
    field = this->mesh.createFieldFromAttachedData<UInt>(field_id, group_name,
                                                         element_kind);
  if (!field)
    field = this->mesh.createFieldFromAttachedData<Real>(field_id, group_name,
                                                         element_kind);

  if (!field)
    AKANTU_DEBUG_WARNING("No field could be found based on name: " << field_id);
  if (field) {
    DumperIOHelper & dumper = mesh.getGroupDumper(dumper_name, group_name);
    this->addDumpGroupFieldToDumper(field_id, field, dumper);
  }

#endif
}

/* -------------------------------------------------------------------------- */

void Model::dump() { mesh.dump(); }

/* -------------------------------------------------------------------------- */

void Model::setDirectory(const std::string & directory) {
  mesh.setDirectory(directory);
}

/* -------------------------------------------------------------------------- */

void Model::setDirectoryToDumper(const std::string & dumper_name,
                                 const std::string & directory) {
  mesh.setDirectoryToDumper(dumper_name, directory);
}

/* -------------------------------------------------------------------------- */

void Model::setTextModeToDumper() { mesh.setTextModeToDumper(); }

/* -------------------------------------------------------------------------- */

} // akantu
