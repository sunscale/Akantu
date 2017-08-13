/**
 * @file   solid_mechanics_model_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Wed Nov 18 2015
 *
 * @brief  Implementation of the inline functions of the SolidMechanicsModel
 * class
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
#include "aka_named_argument.hh"
#include "material_selector.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLID_MECHANICS_MODEL_INLINE_IMPL_CC__
#define __AKANTU_SOLID_MECHANICS_MODEL_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline SolidMechanicsModelOptions::SolidMechanicsModelOptions(
    AnalysisMethod analysis_method)
    : analysis_method(analysis_method) {}

/* -------------------------------------------------------------------------- */
template <typename... pack>
SolidMechanicsModelOptions::SolidMechanicsModelOptions(use_named_args_t,
                                                       pack &&... _pack)
    : SolidMechanicsModelOptions(
          OPTIONAL_NAMED_ARG(analysis_method, _explicit_lumped_mass)) {}

/* -------------------------------------------------------------------------- */
inline Material & SolidMechanicsModel::getMaterial(UInt mat_index) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(mat_index < materials.size(),
                      "The model " << id << " has no material no "
                                   << mat_index);
  AKANTU_DEBUG_OUT();
  return *materials[mat_index];
}

/* -------------------------------------------------------------------------- */
inline const Material & SolidMechanicsModel::getMaterial(UInt mat_index) const {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(mat_index < materials.size(),
                      "The model " << id << " has no material no "
                                   << mat_index);
  AKANTU_DEBUG_OUT();
  return *materials[mat_index];
}

/* -------------------------------------------------------------------------- */
inline Material & SolidMechanicsModel::getMaterial(const std::string & name) {
  AKANTU_DEBUG_IN();
  std::map<std::string, UInt>::const_iterator it =
      materials_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != materials_names_to_id.end(),
                      "The model " << id << " has no material named " << name);
  AKANTU_DEBUG_OUT();
  return *materials[it->second];
}

/* -------------------------------------------------------------------------- */
inline UInt
SolidMechanicsModel::getMaterialIndex(const std::string & name) const {
  AKANTU_DEBUG_IN();
  std::map<std::string, UInt>::const_iterator it =
      materials_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != materials_names_to_id.end(),
                      "The model " << id << " has no material named " << name);
  AKANTU_DEBUG_OUT();
  return it->second;
}

/* -------------------------------------------------------------------------- */
inline const Material &
SolidMechanicsModel::getMaterial(const std::string & name) const {
  AKANTU_DEBUG_IN();
  std::map<std::string, UInt>::const_iterator it =
      materials_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != materials_names_to_id.end(),
                      "The model " << id << " has no material named " << name);
  AKANTU_DEBUG_OUT();
  return *materials[it->second];
}

/* -------------------------------------------------------------------------- */
inline void
SolidMechanicsModel::setMaterialSelector(MaterialSelector & selector) {
  if (is_default_material_selector)
    delete material_selector;
  material_selector = &selector;
  is_default_material_selector = false;
}


/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::splitElementByMaterial(
    const Array<Element> & elements, Array<Element> * elements_per_mat) const {
  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  const Array<UInt> * mat_indexes = NULL;
  const Array<UInt> * mat_loc_num = NULL;

  Array<Element>::const_iterator<Element> it = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    Element el = *it;

    if (el.type != current_element_type ||
        el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type = el.ghost_type;
      mat_indexes = &(this->material_index(el.type, el.ghost_type));
      mat_loc_num = &(this->material_local_numbering(el.type, el.ghost_type));
    }

    UInt old_id = el.element;
    el.element = (*mat_loc_num)(old_id);
    elements_per_mat[(*mat_indexes)(old_id)].push_back(el);
  }
}

/* -------------------------------------------------------------------------- */
inline UInt
SolidMechanicsModel::getNbData(const Array<Element> & elements,
                               const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = 0;

  Array<Element>::const_iterator<Element> it = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

  switch (tag) {
  case _gst_material_id: {
    size += elements.size() * sizeof(UInt);
    break;
  }
  case _gst_smm_mass: {
    size += nb_nodes_per_element * sizeof(Real) *
            Model::spatial_dimension; // mass vector
    break;
  }
  case _gst_smm_for_gradu: {
    size += nb_nodes_per_element * Model::spatial_dimension *
            sizeof(Real); // displacement
    break;
  }
  case _gst_smm_boundary: {
    // force, displacement, boundary
    size += nb_nodes_per_element * Model::spatial_dimension *
            (2 * sizeof(Real) + sizeof(bool));
    break;
  }
  case _gst_for_dump: {
    // displacement, velocity, acceleration, residual, force
    size += nb_nodes_per_element * Model::spatial_dimension * sizeof(Real) * 5;
    break;
  }
  default: {}
  }

  if (tag != _gst_material_id) {
    Array<Element> * elements_per_mat = new Array<Element>[materials.size()];
    this->splitElementByMaterial(elements, elements_per_mat);

    for (UInt i = 0; i < materials.size(); ++i) {
      size += materials[i]->getNbData(elements_per_mat[i], tag);
    }
    delete[] elements_per_mat;
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void
SolidMechanicsModel::packData(CommunicationBuffer & buffer,
                              const Array<Element> & elements,
                              const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  switch (tag) {
  case _gst_material_id: {
    this->packElementalDataHelper(material_index, buffer, elements, false,
                                  getFEEngine());
    break;
  }
  case _gst_smm_mass: {
    packNodalDataHelper(*mass, buffer, elements, mesh);
    break;
  }
  case _gst_smm_for_gradu: {
    packNodalDataHelper(*displacement, buffer, elements, mesh);
    break;
  }
  case _gst_for_dump: {
    packNodalDataHelper(*displacement, buffer, elements, mesh);
    packNodalDataHelper(*velocity, buffer, elements, mesh);
    packNodalDataHelper(*acceleration, buffer, elements, mesh);
    packNodalDataHelper(*internal_force, buffer, elements, mesh);
    packNodalDataHelper(*external_force, buffer, elements, mesh);
    break;
  }
  case _gst_smm_boundary: {
    packNodalDataHelper(*external_force, buffer, elements, mesh);
    packNodalDataHelper(*velocity, buffer, elements, mesh);
    packNodalDataHelper(*blocked_dofs, buffer, elements, mesh);
    break;
  }
  default: {}
  }

  if (tag != _gst_material_id) {
    Array<Element> * elements_per_mat = new Array<Element>[materials.size()];
    splitElementByMaterial(elements, elements_per_mat);

    for (UInt i = 0; i < materials.size(); ++i) {
      materials[i]->packData(buffer, elements_per_mat[i], tag);
    }

    delete[] elements_per_mat;
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::unpackData(CommunicationBuffer & buffer,
                                            const Array<Element> & elements,
                                            const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  switch (tag) {
  case _gst_material_id: {
    unpackElementalDataHelper(material_index, buffer, elements, false,
                              getFEEngine());
    break;
  }
  case _gst_smm_mass: {
    unpackNodalDataHelper(*mass, buffer, elements, mesh);
    break;
  }
  case _gst_smm_for_gradu: {
    unpackNodalDataHelper(*displacement, buffer, elements, mesh);
    break;
  }
  case _gst_for_dump: {
    unpackNodalDataHelper(*displacement, buffer, elements, mesh);
    unpackNodalDataHelper(*velocity, buffer, elements, mesh);
    unpackNodalDataHelper(*acceleration, buffer, elements, mesh);
    unpackNodalDataHelper(*internal_force, buffer, elements, mesh);
    unpackNodalDataHelper(*external_force, buffer, elements, mesh);
    break;
  }
  case _gst_smm_boundary: {
    unpackNodalDataHelper(*external_force, buffer, elements, mesh);
    unpackNodalDataHelper(*velocity, buffer, elements, mesh);
    unpackNodalDataHelper(*blocked_dofs, buffer, elements, mesh);
    break;
  }
  default: {}
  }

  if (tag != _gst_material_id) {
    Array<Element> * elements_per_mat = new Array<Element>[materials.size()];
    splitElementByMaterial(elements, elements_per_mat);

    for (UInt i = 0; i < materials.size(); ++i) {
      materials[i]->unpackData(buffer, elements_per_mat[i], tag);
    }

    delete[] elements_per_mat;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline UInt
SolidMechanicsModel::getNbData(const Array<UInt> & dofs,
                               const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  //  UInt nb_nodes = mesh.getNbNodes();

  switch (tag) {
  case _gst_smm_uv: {
    size += sizeof(Real) * Model::spatial_dimension * 2;
    break;
  }
  case _gst_smm_res: {
    size += sizeof(Real) * Model::spatial_dimension;
    break;
  }
  case _gst_smm_mass: {
    size += sizeof(Real) * Model::spatial_dimension;
    break;
  }
  case _gst_for_dump: {
    size += sizeof(Real) * Model::spatial_dimension * 5;
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size * dofs.size();
}

/* -------------------------------------------------------------------------- */
inline void
SolidMechanicsModel::packData(CommunicationBuffer & buffer,
                              const Array<UInt> & dofs,
                              const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  switch (tag) {
  case _gst_smm_uv: {
    packDOFDataHelper(*displacement, buffer, dofs);
    packDOFDataHelper(*velocity, buffer, dofs);
    break;
  }
  case _gst_smm_res: {
    packDOFDataHelper(*internal_force, buffer, dofs);
    break;
  }
  case _gst_smm_mass: {
    packDOFDataHelper(*mass, buffer, dofs);
    break;
  }
  case _gst_for_dump: {
    packDOFDataHelper(*displacement, buffer, dofs);
    packDOFDataHelper(*velocity, buffer, dofs);
    packDOFDataHelper(*acceleration, buffer, dofs);
    packDOFDataHelper(*internal_force, buffer, dofs);
    packDOFDataHelper(*external_force, buffer, dofs);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::unpackData(CommunicationBuffer & buffer,
                                            const Array<UInt> & dofs,
                                            const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  switch (tag) {
  case _gst_smm_uv: {
    unpackDOFDataHelper(*displacement, buffer, dofs);
    unpackDOFDataHelper(*velocity, buffer, dofs);
    break;
  }
  case _gst_smm_res: {
    unpackDOFDataHelper(*internal_force, buffer, dofs);
    break;
  }
  case _gst_smm_mass: {
    unpackDOFDataHelper(*mass, buffer, dofs);
    break;
  }
  case _gst_for_dump: {
    unpackDOFDataHelper(*displacement, buffer, dofs);
    unpackDOFDataHelper(*velocity, buffer, dofs);
    unpackDOFDataHelper(*acceleration, buffer, dofs);
    unpackDOFDataHelper(*internal_force, buffer, dofs);
    unpackDOFDataHelper(*external_force, buffer, dofs);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
} // namespace akantu

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_INLINE_IMPL_CC__ */
