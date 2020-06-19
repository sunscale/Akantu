/**
 * @file   material_selector_cohesive.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Dec 11 2015
 * @date last modification: Mon Dec 18 2017
 *
 * @brief  Material selector for cohesive elements
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
#include "material_selector_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
DefaultMaterialCohesiveSelector::DefaultMaterialCohesiveSelector(
    const SolidMechanicsModelCohesive & model)
    : facet_material(model.getFacetMaterial()), mesh(model.getMesh()) {
  // backward compatibility v3: to get the former behavior back when the user
  // creates its own selector
  this->fallback_selector =
      std::make_shared<DefaultMaterialSelector>(model.getMaterialByElement());
}

/* -------------------------------------------------------------------------- */
UInt DefaultMaterialCohesiveSelector::operator()(const Element & element) {
  if (Mesh::getKind(element.type) == _ek_cohesive) {
    try {
      const Array<Element> & cohesive_el_to_facet =
          mesh.getMeshFacets().getSubelementToElement(element.type,
                                                      element.ghost_type);
      bool third_dimension = (mesh.getSpatialDimension() == 3);
      const Element & facet =
          cohesive_el_to_facet(element.element, third_dimension);
      if (facet_material.exists(facet.type, facet.ghost_type)) {
        return facet_material(facet.type, facet.ghost_type)(facet.element);
      } else {
        return fallback_value;
      }
    } catch (...) {
      return fallback_value;
    }
  } else if (Mesh::getSpatialDimension(element.type) ==
             mesh.getSpatialDimension() - 1) {
    return facet_material(element.type, element.ghost_type)(element.element);
  } else {
    return MaterialSelector::operator()(element);
  }
}

/* -------------------------------------------------------------------------- */
MeshDataMaterialCohesiveSelector::MeshDataMaterialCohesiveSelector(
    const SolidMechanicsModelCohesive & model)
    : model(model), mesh_facets(model.getMeshFacets()),
      material_index(mesh_facets.getData<std::string>("physical_names")) {
  third_dimension = (model.getSpatialDimension() == 3);
  // backward compatibility v3: to get the former behavior back when the user
  // creates its own selector
  this->fallback_selector =
      std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names",
                                                              model);
}

/* -------------------------------------------------------------------------- */
UInt MeshDataMaterialCohesiveSelector::operator()(const Element & element) {
  if (Mesh::getKind(element.type) == _ek_cohesive or
      Mesh::getSpatialDimension(element.type) == mesh_facets.getSpatialDimension() - 1) {
    Element facet;
    if (Mesh::getKind(element.type) == _ek_cohesive) {
      facet = mesh_facets.getSubelementToElement(element.type,
						 element.ghost_type)(element.element,
								     third_dimension);
    } else {
      facet = element;
    }
  
    try {
      std::string material_name = this->material_index(facet);
      return this->model.getMaterialIndex(material_name);
    } catch (...) {
      return fallback_value;
    }
  }
  return MaterialSelector::operator()(element);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
MaterialCohesiveRulesSelector::MaterialCohesiveRulesSelector(
    const SolidMechanicsModelCohesive & model,
    const MaterialCohesiveRules & rules,
    ID mesh_data_id) // what we have here is the name of model and also
                     // the name of different materials
    : model(model), mesh_data_id(std::move(mesh_data_id)),
      mesh(model.getMesh()), mesh_facets(model.getMeshFacets()),
      spatial_dimension(model.getSpatialDimension()), rules(rules) {

  // cohesive fallback
  this->fallback_selector =
      std::make_shared<DefaultMaterialCohesiveSelector>(model);

  // non cohesive fallback
  this->fallback_selector->setFallback(
      std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names",
                                                              model));
}

/* -------------------------------------------------------------------------- */
UInt MaterialCohesiveRulesSelector::operator()(const Element & element) {
  if (mesh_facets.getSpatialDimension(element.type) ==
      (spatial_dimension - 1)) {
    const std::vector<Element> & element_to_subelement =
        mesh_facets.getElementToSubelement(element.type,
                                           element.ghost_type)(element.element);
    // Array<bool> & facets_check = model.getFacetsCheck();

    const Element & el1 = element_to_subelement[0];
    const Element & el2 = element_to_subelement[1];

    ID id1 = mesh.getData<std::string>(mesh_data_id, el1.type,
                                       el1.ghost_type)(el1.element);

    ID id2 = id1;
    if (el2 != ElementNull) {
      id2 = mesh.getData<std::string>(mesh_data_id, el2.type,
                                      el2.ghost_type)(el2.element);
    }

    auto rit = rules.find(std::make_pair(id1, id2));
    if (rit == rules.end()) {
      rit = rules.find(std::make_pair(id2, id1));
    }

    if (rit != rules.end()) {
      return model.getMaterialIndex(rit->second);
    }
  }

  return MaterialSelector::operator()(element);
}

} // namespace akantu
