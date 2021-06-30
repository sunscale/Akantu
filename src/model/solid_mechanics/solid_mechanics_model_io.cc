/**
 * @file   solid_mechanics_model_io.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Jul 09 2017
 * @date last modification: Sun Dec 03 2017
 *
 * @brief  Dumpable part of the SolidMechnicsModel
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"

#include "group_manager_inline_impl.hh"

#include "dumpable_inline_impl.hh"
#ifdef AKANTU_USE_IOHELPER
#include "dumper_element_partition.hh"
#include "dumper_elemental_field.hh"
#include "dumper_field.hh"
#include "dumper_homogenizing_field.hh"
#include "dumper_internal_material_field.hh"
#include "dumper_iohelper.hh"
#include "dumper_material_padders.hh"
#include "dumper_paraview.hh"
#endif

namespace akantu {

/* -------------------------------------------------------------------------- */
bool SolidMechanicsModel::isInternal(const std::string & field_name,
                                     ElementKind element_kind) {
  /// check if at least one material contains field_id as an internal
  for (auto & material : materials) {
    bool is_internal = material->isInternal<Real>(field_name, element_kind);
    if (is_internal) {
      return true;
    }
  }

  return false;
}

/* -------------------------------------------------------------------------- */
ElementTypeMap<UInt>
SolidMechanicsModel::getInternalDataPerElem(const std::string & field_name,
                                            ElementKind element_kind) {

  if (!(this->isInternal(field_name, element_kind))) {
    AKANTU_EXCEPTION("unknown internal " << field_name);
  }

  for (auto & material : materials) {
    if (material->isInternal<Real>(field_name, element_kind)) {
      return material->getInternalDataPerElem<Real>(field_name, element_kind);
    }
  }

  return ElementTypeMap<UInt>();
}

/* -------------------------------------------------------------------------- */
ElementTypeMapArray<Real> &
SolidMechanicsModel::flattenInternal(const std::string & field_name,
                                     ElementKind kind,
                                     const GhostType ghost_type) {
  auto key = std::make_pair(field_name, kind);

  ElementTypeMapArray<Real> * internal_flat;

  auto it = this->registered_internals.find(key);
  if (it == this->registered_internals.end()) {
    auto internal = std::make_unique<ElementTypeMapArray<Real>>(
        field_name, this->id);

    internal_flat = internal.get();
    this->registered_internals[key] = std::move(internal);
  } else {
    internal_flat = it->second.get();
  }

  for (auto type :
       mesh.elementTypes(Model::spatial_dimension, ghost_type, kind)) {
    if (internal_flat->exists(type, ghost_type)) {
      auto & internal = (*internal_flat)(type, ghost_type);
      internal.resize(0);
    }
  }

  for (auto & material : materials) {
    if (material->isInternal<Real>(field_name, kind)) {
      material->flattenInternal(field_name, *internal_flat, ghost_type, kind);
    }
  }

  return *internal_flat;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::flattenAllRegisteredInternals(
    ElementKind kind) {
  ElementKind _kind;
  ID _id;

  for (auto & internal : this->registered_internals) {
    std::tie(_id, _kind) = internal.first;
    if (kind == _kind) {
      this->flattenInternal(_id, kind);
    }
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onDump() {
  this->flattenAllRegisteredInternals(_ek_regular);
}

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
std::shared_ptr<dumpers::Field> SolidMechanicsModel::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag, UInt spatial_dimension,
    ElementKind kind) {

  std::shared_ptr<dumpers::Field> field;

  if (field_name == "partitions") {
    field = mesh.createElementalField<UInt, dumpers::ElementPartitionField>(
        mesh.getConnectivities(), group_name, spatial_dimension, kind);
  } else if (field_name == "material_index") {
    field = mesh.createElementalField<UInt, Vector, dumpers::ElementalField>(
        material_index, group_name, spatial_dimension, kind);
  } else {
    // this copy of field_name is used to compute derivated data such as
    // strain and von mises stress that are based on grad_u and stress
    std::string field_name_copy(field_name);

    if (field_name == "strain" || field_name == "Green strain" ||
        field_name == "principal strain" ||
        field_name == "principal Green strain") {
      field_name_copy = "grad_u";
    } else if (field_name == "Von Mises stress") {
      field_name_copy = "stress";
    }

    bool is_internal = this->isInternal(field_name_copy, kind);

    if (is_internal) {
      auto nb_data_per_elem =
          this->getInternalDataPerElem(field_name_copy, kind);
      auto & internal_flat = this->flattenInternal(field_name_copy, kind);

      field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
          internal_flat, group_name, spatial_dimension, kind, nb_data_per_elem);

      std::unique_ptr<dumpers::ComputeFunctorInterface> func;
      if (field_name == "strain") {
        func = std::make_unique<dumpers::ComputeStrain<false>>(*this);
      } else if (field_name == "Von Mises stress") {
        func = std::make_unique<dumpers::ComputeVonMisesStress>(*this);
      } else if (field_name == "Green strain") {
        func = std::make_unique<dumpers::ComputeStrain<true>>(*this);
      } else if (field_name == "principal strain") {
        func = std::make_unique<dumpers::ComputePrincipalStrain<false>>(*this);
      } else if (field_name == "principal Green strain") {
        func = std::make_unique<dumpers::ComputePrincipalStrain<true>>(*this);
      }

      if (func) {
        field = dumpers::FieldComputeProxy::createFieldCompute(field,
                                                              std::move(func));
      }
      // treat the paddings
      if (padding_flag) {
        if (field_name == "stress") {
          if (spatial_dimension == 2) {
            auto foo = std::make_unique<dumpers::StressPadder<2>>(*this);
            field = dumpers::FieldComputeProxy::createFieldCompute(
                field, std::move(foo));
          }
        } else if (field_name == "strain" || field_name == "Green strain") {
          if (spatial_dimension == 2) {
            auto foo = std::make_unique<dumpers::StrainPadder<2>>(*this);
            field = dumpers::FieldComputeProxy::createFieldCompute(
                field, std::move(foo));
          }
        }
      }

      // homogenize the field
      auto foo = dumpers::HomogenizerProxy::createHomogenizer(*field);

      field =
          dumpers::FieldComputeProxy::createFieldCompute(field, std::move(foo));
    }
  }
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field>
SolidMechanicsModel::createNodalFieldReal(const std::string & field_name,
                                          const std::string & group_name,
                                          bool padding_flag) {

  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["displacement"] = this->displacement.get();
  real_nodal_fields["mass"] = this->mass.get();
  real_nodal_fields["velocity"] = this->velocity.get();
  real_nodal_fields["acceleration"] = this->acceleration.get();
  real_nodal_fields["external_force"] = this->external_force.get();
  real_nodal_fields["internal_force"] = this->internal_force.get();
  real_nodal_fields["increment"] = this->displacement_increment.get();

  if (field_name == "force") {
    AKANTU_EXCEPTION("The 'force' field has been renamed in 'external_force'");
  } else if (field_name == "residual") {
    AKANTU_EXCEPTION(
        "The 'residual' field has been replaced by 'internal_force'");
  }

  std::shared_ptr<dumpers::Field> field;
  if (padding_flag) {
    field = this->mesh.createNodalField(real_nodal_fields[field_name],
                                        group_name, 3);
  } else {
    field =
        this->mesh.createNodalField(real_nodal_fields[field_name], group_name);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> SolidMechanicsModel::createNodalFieldBool(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  std::map<std::string, Array<bool> *> uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"] = blocked_dofs.get();

  std::shared_ptr<dumpers::Field> field;
  field = mesh.createNodalField(uint_nodal_fields[field_name], group_name);
  return field;
}
/* -------------------------------------------------------------------------- */
#else
/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field>
SolidMechanicsModel::createElementalField(const std::string &,
                                          const std::string &, bool,
                                          const UInt &, ElementKind) {
  return nullptr;
}
/* --------------------------------------------------------------------------
 */
std::shaed_ptr<dumpers::Field>
SolidMechanicsModel::createNodalFieldReal(const std::string &,
                                          const std::string &, bool) {
  return nullptr;
}

/* --------------------------------------------------------------------------
 */
std::shared_ptr<dumpers::Field>
SolidMechanicsModel::createNodalFieldBool(const std::string &,
                                          const std::string &, bool) {
  return nullptr;
}

#endif
/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::dump(const std::string & dumper_name) {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  mesh.dump(dumper_name);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::dump(const std::string & dumper_name, UInt step) {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  mesh.dump(dumper_name, step);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::dump(const std::string & dumper_name, Real time,
                               UInt step) {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  mesh.dump(dumper_name, time, step);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::dump() {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  mesh.dump();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::dump(UInt step) {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  mesh.dump(step);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::dump(Real time, UInt step) {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  mesh.dump(time, step);
}

} // namespace akantu
