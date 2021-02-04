/**
 * @file   model.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Oct 03 2011
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  implementation of model common parts
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "model.hh"
#include "communicator.hh"
#include "data_accessor.hh"
#include "element_group.hh"
#include "element_synchronizer.hh"
#include "synchronizer_registry.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Model::Model(Mesh & mesh, const ModelType & type, UInt dim, const ID & id,
             const MemoryID & memory_id)
    : Memory(id, memory_id), ModelSolver(mesh, type, id, memory_id), mesh(mesh),
      spatial_dimension(dim == _all_dimensions ? mesh.getSpatialDimension()
                                               : dim),
      parser(getStaticParser()) {
  AKANTU_DEBUG_IN();

  this->mesh.registerEventHandler(*this, _ehp_model);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Model::~Model() = default;

/* -------------------------------------------------------------------------- */
void Model::initFullImpl(const ModelOptions & options) {
  AKANTU_DEBUG_IN();

  method = options.analysis_method;
  if (!this->hasDefaultSolver()) {
    this->initNewSolver(this->method);
  }
  initModel();

  initFEEngineBoundary();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Model::initNewSolver(const AnalysisMethod & method) {
  ID solver_name;
  TimeStepSolverType tss_type;
  std::tie(solver_name, tss_type) = this->getDefaultSolverID(method);

  if (not this->hasSolver(solver_name)) {
    ModelSolverOptions options = this->getDefaultSolverOptions(tss_type);
    this->getNewSolver(solver_name, tss_type, options.non_linear_solver_type);

    for (auto && is_type : options.integration_scheme_type) {
      if (!this->hasIntegrationScheme(solver_name, is_type.first)) {
        this->setIntegrationScheme(solver_name, is_type.first, is_type.second,
                                   options.solution_type[is_type.first]);
      }
    }
  }

  this->method = method;
  this->setDefaultSolver(solver_name);
}

/* -------------------------------------------------------------------------- */
void Model::initFEEngineBoundary() {
  FEEngine & fem_boundary = getFEEngineBoundary();
  fem_boundary.initShapeFunctions(_not_ghost);
  fem_boundary.initShapeFunctions(_ghost);

  fem_boundary.computeNormalsOnIntegrationPoints(_not_ghost);
  fem_boundary.computeNormalsOnIntegrationPoints(_ghost);
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
  for (auto & group : mesh.iterateElementGroups()) {
    group.dump();
  }
}

/* -------------------------------------------------------------------------- */
void Model::setGroupDirectory(const std::string & directory) {
  for (auto & group : mesh.iterateElementGroups()) {
    group.setDirectory(directory);
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
                                      std::shared_ptr<dumpers::Field> field,
                                      DumperIOHelper & dumper) {
#ifdef AKANTU_USE_IOHELPER
  dumper.registerField(field_id, std::move(field));
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
  this->addDumpGroupFieldToDumper(dumper_name, field_id, "all",
                                  dumper_default_element_kind, false);
}

/* -------------------------------------------------------------------------- */
void Model::addDumpGroupField(const std::string & field_id,
                              const std::string & group_name) {
  ElementGroup & group = mesh.getElementGroup(group_name);
  this->addDumpGroupFieldToDumper(group.getDefaultDumperName(), field_id,
                                  group_name, dumper_default_element_kind,
                                  false);
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
  this->addDumpGroupFieldToDumper(dumper_name, field_id, "all",
                                  dumper_default_element_kind, true);
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
                                  dumper_default_element_kind, true);
}
/* -------------------------------------------------------------------------- */

void Model::addDumpFieldTensorToDumper(const std::string & dumper_name,
                                       const std::string & field_id) {
  this->addDumpGroupFieldToDumper(dumper_name, field_id, "all",
                                  dumper_default_element_kind, true);
}

/* -------------------------------------------------------------------------- */
void Model::addDumpGroupFieldToDumper(const std::string & dumper_name,
                                      const std::string & field_id,
                                      const std::string & group_name,
                                      ElementKind element_kind,
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
                                      ElementKind element_kind,
                                      bool padding_flag) {

#ifdef AKANTU_USE_IOHELPER
  std::shared_ptr<dumpers::Field> field;

  if (!field) {
    field = this->createNodalFieldReal(field_id, group_name, padding_flag);
  }
  if (!field) {
    field = this->createNodalFieldUInt(field_id, group_name, padding_flag);
  }
  if (!field) {
    field = this->createNodalFieldBool(field_id, group_name, padding_flag);
  }
  if (!field) {
    field = this->createElementalField(field_id, group_name, padding_flag,
                                       spatial_dimension, element_kind);
  }
  if (!field) {
    field = this->mesh.createFieldFromAttachedData<UInt>(field_id, group_name,
                                                         element_kind);
  }
  if (!field) {
    field = this->mesh.createFieldFromAttachedData<Real>(field_id, group_name,
                                                         element_kind);
  }

#ifndef AKANTU_NDEBUG
  if (!field) {
    AKANTU_DEBUG_WARNING("No field could be found based on name: " << field_id);
  }
#endif
  if (field) {
    DumperIOHelper & dumper = mesh.getGroupDumper(dumper_name, group_name);
    this->addDumpGroupFieldToDumper(field_id, field, dumper);
  }

#endif
}

/* -------------------------------------------------------------------------- */
void Model::dump() { mesh.dump(); }

/* -------------------------------------------------------------------------- */
void Model::dump(UInt step) { mesh.dump(step); }

/* -------------------------------------------------------------------------- */
void Model::dump(Real time, UInt step) { mesh.dump(time, step); }

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

} // namespace akantu
