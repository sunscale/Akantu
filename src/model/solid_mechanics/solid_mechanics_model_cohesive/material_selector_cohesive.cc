/**
 * @file   material_selector_cohesive.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Dec 11 2015
 * @date last modification: Mon Dec 14 2015
 *
 * @brief  Material selector for cohesive elements
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "material_selector_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
DefaultMaterialCohesiveSelector::DefaultMaterialCohesiveSelector(
    const SolidMechanicsModelCohesive & model)
    : facet_material(model.getFacetMaterial()), mesh(model.getMesh()) {}

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
    : MeshDataMaterialSelector<std::string>("physical_names", model),
      mesh_facets(model.getMeshFacets()),
      material_index(mesh_facets.getData<UInt>("physical_names")) {
  third_dimension = (model.getSpatialDimension() == 3);
}

/* -------------------------------------------------------------------------- */
UInt MeshDataMaterialCohesiveSelector::operator()(const Element & element) {
  if (Mesh::getKind(element.type) == _ek_cohesive) {
    const Array<Element> & cohesive_el_to_facet =
        mesh_facets.getSubelementToElement(element.type, element.ghost_type);
    const Element & facet =
        cohesive_el_to_facet(element.element, third_dimension);
    UInt material_id =
        material_index(facet.type, facet.ghost_type)(facet.element);
    UInt fallback_id = MaterialSelector::operator()(element);

    return std::max(material_id, fallback_id);
  } else
    return MeshDataMaterialSelector<std::string>::operator()(element);
}

} // namespace akantu
