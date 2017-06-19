/**
 * @file   material_selector_cohesive.hh
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Dec 11 2015
 * @date last modification: Mon Dec 14 2015
 *
 * @brief  Material selectors for cohesive elements
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
#include "material_selector.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
  class SolidMechanicsModelCohesive;
}

namespace akantu {

#ifndef __AKANTU_MATERIAL_SELECTOR_COHESIVE_HH__
#define __AKANTU_MATERIAL_SELECTOR_COHESIVE_HH__

/* -------------------------------------------------------------------------- */
/**
 * class that assigns the first cohesive material by default to the
 * cohesive elements
 */
class DefaultMaterialCohesiveSelector : public DefaultMaterialSelector {
public:
  DefaultMaterialCohesiveSelector(const SolidMechanicsModelCohesive & model);
  virtual UInt operator()(const Element & element);

private:
  const ElementTypeMapArray<UInt> & facet_material;
  const Mesh & mesh;
};

/* -------------------------------------------------------------------------- */
/// To be used with intrinsic elements inserted along mesh physical surfaces
class MeshDataMaterialCohesiveSelector
    : public MeshDataMaterialSelector<std::string> {
public:
  MeshDataMaterialCohesiveSelector(const SolidMechanicsModelCohesive & model);
  virtual UInt operator()(const Element & element);
protected:
  const Mesh &mesh_facets;
  const ElementTypeMapArray<UInt> & material_index;
  bool third_dimension;
};

#endif /* __AKANTU_MATERIAL_SELECTOR_COHESIVE_HH__ */

} // akantu
