/**
 * @file   solid_mechanics_model_material.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Nov 26 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  instatiation of materials
 *
 * @section LICENSE
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
#include "aka_factory.hh"
#include "aka_math.hh"
#include "material_non_local.hh"
#include "mesh_iterators.hh"
#include "non_local_manager.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Material &
SolidMechanicsModel::registerNewMaterial(const ParserSection & section) {
  std::string mat_name;
  std::string mat_type = section.getName();
  std::string opt_param = section.getOption();

  try {
    std::string tmp = section.getParameter("name");
    mat_name = tmp; /** this can seam weird, but there is an ambiguous operator
                     * overload that i couldn't solve. @todo remove the
                     * weirdness of this code
                     */
  } catch (debug::Exception &) {
    AKANTU_ERROR("A material of type \'"
                 << mat_type
                 << "\' in the input file has been defined without a name!");
  }
  Material & mat = this->registerNewMaterial(mat_name, mat_type, opt_param);

  mat.parseSection(section);

  return mat;
}

/* -------------------------------------------------------------------------- */
Material & SolidMechanicsModel::registerNewMaterial(const ID & mat_name,
                                                    const ID & mat_type,
                                                    const ID & opt_param) {
  AKANTU_DEBUG_ASSERT(materials_names_to_id.find(mat_name) ==
                          materials_names_to_id.end(),
                      "A material with this name '"
                          << mat_name << "' has already been registered. "
                          << "Please use unique names for materials");

  UInt mat_count = materials.size();
  materials_names_to_id[mat_name] = mat_count;

  std::stringstream sstr_mat;
  sstr_mat << this->id << ":" << mat_count << ":" << mat_type;
  ID mat_id = sstr_mat.str();

  std::unique_ptr<Material> material = MaterialFactory::getInstance().allocate(
      mat_type, spatial_dimension, opt_param, *this, mat_id);

  materials.push_back(std::move(material));

  return *(materials.back());
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::instantiateMaterials() {
  ParserSection model_section;
  bool is_empty;
  std::tie(model_section, is_empty) = this->getParserSection();

  if (not is_empty) {
    auto model_materials = model_section.getSubSections(ParserType::_material);
    for (const auto & section : model_materials) {
      this->registerNewMaterial(section);
    }
  }

  auto sub_sections = this->parser.getSubSections(ParserType::_material);
  for (const auto & section : sub_sections) {
    this->registerNewMaterial(section);
  }

#ifdef AKANTU_DAMAGE_NON_LOCAL
  for (auto & material : materials) {
    if (dynamic_cast<MaterialNonLocalInterface *>(material.get()) == nullptr)
      continue;

    this->non_local_manager = std::make_unique<NonLocalManager>(
        *this, *this, id + ":non_local_manager", memory_id);
    break;
  }
#endif

  if (materials.empty())
    AKANTU_EXCEPTION("No materials where instantiated for the model"
                     << getID());
  are_materials_instantiated = true;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assignMaterialToElements(
    const ElementTypeMapArray<UInt> * filter) {

  for_each_element(
      mesh,
      [&](auto && element) {
        UInt mat_index = (*material_selector)(element);
        AKANTU_DEBUG_ASSERT(
            mat_index < materials.size(),
            "The material selector returned an index that does not exists");
        material_index(element) = mat_index;
      },
      _element_filter = filter, _ghost_type = _not_ghost);

  if (non_local_manager)
    non_local_manager->synchronize(*this, SynchronizationTag::_material_id);

  for_each_element(mesh,
                   [&](auto && element) {
                     auto mat_index = material_index(element);
                     auto index = materials[mat_index]->addElement(element);
                     material_local_numbering(element) = index;
                   },
                   _element_filter = filter, _ghost_type = _not_ghost);

  // synchronize the element material arrays
  this->synchronize(SynchronizationTag::_material_id);
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

  this->synchronize(SynchronizationTag::_smm_init_mat);

  if (this->non_local_manager) {
    this->non_local_manager->initialize();
  }
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

    for (auto type :
         mesh.elementTypes(spatial_dimension, ghost_type, _ek_not_defined)) {
      element.type = type;

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
                          spatial_dimension * spatial_dimension,
                      "The prescribed grad_u is not of the good size");
  for (auto & material : materials) {
    if (material->getName() == material_name)
      material->applyEigenGradU(prescribed_eigen_grad_u, ghost_type);
  }
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
