/**
 * @file   solid_mechanics_model_material.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Nov 26 2010
 * @date last modification: Mon Nov 16 2015
 *
 * @brief  instatiation of materials
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
#include "aka_math.hh"
#include "material_list.hh"
#include "solid_mechanics_model.hh"
#ifdef AKANTU_DAMAGE_NON_LOCAL
#include "non_local_manager.hh"
#endif
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
#define AKANTU_INTANTIATE_MATERIAL_BY_DIM_NO_TMPL(dim, elem)                   \
  registerNewMaterial<BOOST_PP_ARRAY_ELEM(1, elem) < dim>> (section)

#define AKANTU_INTANTIATE_MATERIAL_BY_DIM_TMPL_EACH(r, data, i, elem)          \
  BOOST_PP_EXPR_IF(BOOST_PP_NOT_EQUAL(0, i), else)                             \
  if (opt_param == BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2, 0, elem))) {      \
    registerNewMaterial<BOOST_PP_ARRAY_ELEM(1, data) <                         \
                            BOOST_PP_ARRAY_ELEM(0, data),                      \
                        BOOST_PP_SEQ_ENUM(BOOST_PP_TUPLE_ELEM(2, 1, elem))>>   \
        (section);                                                             \
  }

#define AKANTU_INTANTIATE_MATERIAL_BY_DIM_TMPL(dim, elem)                      \
  BOOST_PP_SEQ_FOR_EACH_I(AKANTU_INTANTIATE_MATERIAL_BY_DIM_TMPL_EACH,         \
                          (2, (dim, BOOST_PP_ARRAY_ELEM(1, elem))),            \
                          BOOST_PP_ARRAY_ELEM(2, elem))                        \
  else {                                                                       \
    AKANTU_INTANTIATE_MATERIAL_BY_DIM_NO_TMPL(dim, elem);                      \
  }

#define AKANTU_INTANTIATE_MATERIAL_BY_DIM(dim, elem)                           \
  BOOST_PP_IF(BOOST_PP_EQUAL(3, BOOST_PP_ARRAY_SIZE(elem)),                    \
              AKANTU_INTANTIATE_MATERIAL_BY_DIM_TMPL,                          \
              AKANTU_INTANTIATE_MATERIAL_BY_DIM_NO_TMPL)                       \
  (dim, elem)

#define AKANTU_INTANTIATE_MATERIAL(elem)                                       \
  switch (Model::spatial_dimension) {                                          \
  case 1: {                                                                    \
    AKANTU_INTANTIATE_MATERIAL_BY_DIM(1, elem);                                \
    break;                                                                     \
  }                                                                            \
  case 2: {                                                                    \
    AKANTU_INTANTIATE_MATERIAL_BY_DIM(2, elem);                                \
    break;                                                                     \
  }                                                                            \
  case 3: {                                                                    \
    AKANTU_INTANTIATE_MATERIAL_BY_DIM(3, elem);                                \
    break;                                                                     \
  }                                                                            \
  }

#define AKANTU_INTANTIATE_MATERIAL_IF(elem)                                    \
  if (mat_type == BOOST_PP_STRINGIZE(BOOST_PP_ARRAY_ELEM(0, elem))) {          \
    AKANTU_INTANTIATE_MATERIAL(elem);                                          \
  }

#define AKANTU_INTANTIATE_OTHER_MATERIAL(r, data, elem)                        \
  else AKANTU_INTANTIATE_MATERIAL_IF(elem)

#define AKANTU_INSTANTIATE_MATERIALS()                                         \
  do {                                                                         \
    AKANTU_INTANTIATE_MATERIAL_IF(BOOST_PP_SEQ_HEAD(AKANTU_MATERIAL_LIST))     \
    BOOST_PP_SEQ_FOR_EACH(AKANTU_INTANTIATE_OTHER_MATERIAL, _,                 \
                          BOOST_PP_SEQ_TAIL(AKANTU_MATERIAL_LIST))             \
    else {                                                                     \
      if (getStaticParser().isPermissive())                                    \
        AKANTU_DEBUG_INFO("Malformed material file "                           \
                          << ": unknown material type '" << mat_type << "'");  \
      else                                                                     \
        AKANTU_DEBUG_WARNING("Malformed material file "                        \
                             << ": unknown material type " << mat_type         \
                             << ". This is perhaps a user"                     \
                             << " defined material ?");                        \
    }                                                                          \
  } while (0)

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::instantiateMaterials() {
  std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
      sub_sect = this->parser->getSubSections(_st_material);

  Parser::const_section_iterator it = sub_sect.first;
  for (; it != sub_sect.second; ++it) {
    const ParserSection & section = *it;
    std::string mat_type = section.getName();
    std::string opt_param = section.getOption();
    AKANTU_INSTANTIATE_MATERIALS();
  }

  are_materials_instantiated = true;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assignMaterialToElements(
    const ElementTypeMapArray<UInt> * filter) {

  Element element;
  element.ghost_type = _not_ghost;

  auto element_types =
      mesh.elementTypes(Model::spatial_dimension, _not_ghost, _ek_not_defined);
  if (filter != NULL) {
    element_types = filter->elementTypes(Model::spatial_dimension, _not_ghost,
                                         _ek_not_defined);
  }

  // Fill the element material array from the material selector
  for (auto type : element_types) {
    UInt nb_element = mesh.getNbElement(type, _not_ghost);

    const Array<UInt> * filter_array = NULL;
    if (filter != NULL) {
      filter_array = &((*filter)(type, _not_ghost));
      nb_element = filter_array->getSize();
    }

    element.type = type;
    element.kind = mesh.getKind(element.type);
    Array<UInt> & mat_indexes = material_index(type, _not_ghost);
    for (UInt el = 0; el < nb_element; ++el) {
      if (filter != NULL)
        element.element = (*filter_array)(el);
      else
        element.element = el;

      UInt mat_index = (*material_selector)(element);
      AKANTU_DEBUG_ASSERT(
          mat_index < materials.size(),
          "The material selector returned an index that does not exists");
      mat_indexes(element.element) = mat_index;
    }
  }

  // synchronize the element material arrays
  this->synchronize(_gst_material_id);

  /// fill the element filters of the materials using the element_material
  /// arrays
  for (auto ghost_type : ghost_types) {
    element_types = mesh.elementTypes(Model::spatial_dimension, ghost_type,
                                      _ek_not_defined);

    if (filter != NULL) {
      element_types = filter->elementTypes(Model::spatial_dimension, ghost_type,
                                           _ek_not_defined);
    }

    for (auto type : element_types) {
      UInt nb_element = mesh.getNbElement(type, ghost_type);

      const Array<UInt> * filter_array = NULL;
      if (filter != NULL) {
        filter_array = &((*filter)(type, ghost_type));
        nb_element = filter_array->getSize();
      }

      Array<UInt> & mat_indexes = material_index(type, ghost_type);
      Array<UInt> & mat_local_num = material_local_numbering(type, ghost_type);
      for (UInt el = 0; el < nb_element; ++el) {
        UInt element;

        if (filter != NULL)
          element = (*filter_array)(el);
        else
          element = el;

        UInt mat_index = mat_indexes(element);
        UInt index =
            materials[mat_index]->addElement(type, element, ghost_type);
        mat_local_num(element) = index;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initMaterials() {
  AKANTU_DEBUG_ASSERT(materials.size() != 0, "No material to initialize !");

  if (!are_materials_instantiated)
    instantiateMaterials();

  this->assignMaterialToElements();

  for (auto & material : materials) {
    /// init internals properties
    material->initMaterial();
  }

  this->synchronize(_gst_smm_init_mat);

  // initialize mass
  switch (method) {
  case _explicit_lumped_mass:
    assembleMassLumped();
    break;
  case _explicit_consistent_mass:
  case _implicit_dynamic:
    assembleMass();
    break;
  case _static:
    break;
  default:
    AKANTU_EXCEPTION("analysis method not recognised by SolidMechanicsModel");
    break;
  }

// initialize the previous displacement array if at least on material needs it
// for (auto & material : materials) {
//   if (material->isFiniteDeformation() || material->isInelasticDeformation())
//   {
//     initArraysPreviousDisplacment();
//     break;
//   }
// }

#ifdef AKANTU_DAMAGE_NON_LOCAL
  /// initialize the non-local manager for non-local computations
  this->initializeNonLocal();
#endif
}

/* -------------------------------------------------------------------------- */
Int SolidMechanicsModel::getInternalIndexFromID(const ID & id) const {
  AKANTU_DEBUG_IN();

  auto it = materials.begin();
  auto end = materials.end();

  for (; it != end; ++it)
    if ((*it)->getID() == id) {
      AKANTU_DEBUG_OUT();
      return (it - materials.begin());
    }

  AKANTU_DEBUG_OUT();
  return -1;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::reassignMaterial() {
  AKANTU_DEBUG_IN();

  std::vector<Array<Element>> element_to_add(materials.size());
  std::vector<Array<Element>> element_to_remove(materials.size());

  Element element;
  for (auto ghost_type : ghost_types) {
    element.ghost_type = ghost_type;

    for (auto type : mesh.elementTypes(Model::spatial_dimension, ghost_type,
                                       _ek_not_defined)) {
      element.type = type;
      element.kind = Mesh::getKind(type);

      UInt nb_element = mesh.getNbElement(type, ghost_type);
      Array<UInt> & mat_indexes = material_index(type, ghost_type);

      for (UInt el = 0; el < nb_element; ++el) {
        element.element = el;

        UInt old_material = mat_indexes(el);
        UInt new_material = (*material_selector)(element);

        if (old_material != new_material) {
          element_to_add[new_material].push_back(element);
          element_to_remove[old_material].push_back(element);
        }
      }
    }
  }

  UInt mat_index = 0;
  for (auto mat_it = materials.begin(); mat_it != materials.end();
       ++mat_it, ++mat_index) {
    (*mat_it)->removeElements(element_to_remove[mat_index]);
    (*mat_it)->addElements(element_to_add[mat_index]);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::applyEigenGradU(
    const Matrix<Real> & prescribed_eigen_grad_u, const ID & material_name,
    const GhostType ghost_type) {

  AKANTU_DEBUG_ASSERT(prescribed_eigen_grad_u.size() ==
                          Model::spatial_dimension * Model::spatial_dimension,
                      "The prescribed grad_u is not of the good size");
  for (auto & material : materials) {
    if (material->getName() == material_name)
      material->applyEigenGradU(prescribed_eigen_grad_u, ghost_type);
  }
}

/* -------------------------------------------------------------------------- */

} // akantu
