/**
 * @file   material_selector_cohesive.hh
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Dec 11 2015
 * @date last modification: Mon Dec 18 2017
 *
 * @brief  Material selectors for cohesive elements
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
#include "material_selector.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

namespace akantu {
class SolidMechanicsModelCohesive;
}

namespace akantu {

#ifndef AKANTU_MATERIAL_SELECTOR_COHESIVE_HH_
#define AKANTU_MATERIAL_SELECTOR_COHESIVE_HH_

/* -------------------------------------------------------------------------- */
/**
 * class that assigns the first cohesive material by default to the
 * cohesive elements
 */
class DefaultMaterialCohesiveSelector : public MaterialSelector {
public:
  DefaultMaterialCohesiveSelector(const SolidMechanicsModelCohesive & model);
  UInt operator()(const Element & element) override;

private:
  const ElementTypeMapArray<UInt> & facet_material;
  const Mesh & mesh;
};

/* -------------------------------------------------------------------------- */
/// To be used with intrinsic elements inserted along mesh physical surfaces
class MeshDataMaterialCohesiveSelector : public MaterialSelector {
public:
  MeshDataMaterialCohesiveSelector(const SolidMechanicsModelCohesive & model);
  UInt operator()(const Element & element) override;

protected:
  const SolidMechanicsModelCohesive & model;
  const Mesh & mesh_facets;
  const ElementTypeMapArray<std::string> & material_index;
  bool third_dimension;
};

/// bulk1, bulk2 -> cohesive
using MaterialCohesiveRules = std::map<std::pair<ID, ID>, ID>;

/* -------------------------------------------------------------------------- */
class MaterialCohesiveRulesSelector : public MaterialSelector {
public:
  MaterialCohesiveRulesSelector(const SolidMechanicsModelCohesive & model,
                                const MaterialCohesiveRules & rules,
                                ID mesh_data_id = "physical_names");
  UInt operator()(const Element & element) override;

private:
  const SolidMechanicsModelCohesive & model;
  ID mesh_data_id;
  const Mesh & mesh;
  const Mesh & mesh_facets;
  UInt spatial_dimension;
  MaterialCohesiveRules rules;
};

#endif /* AKANTU_MATERIAL_SELECTOR_COHESIVE_HH_ */

} // namespace akantu
