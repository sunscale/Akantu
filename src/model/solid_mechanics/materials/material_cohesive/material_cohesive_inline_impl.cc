/**
 * @file   material_cohesive_inline_impl.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 23 2012
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  MaterialCohesive inline implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
inline UInt MaterialCohesive::addFacet(const Element & element) {
  Array<UInt> & f_filter = facet_filter(element.type, element.ghost_type);
  f_filter.push_back(element.element);
  return f_filter.getSize()-1;
}


/* -------------------------------------------------------------------------- */
template<ElementType type>
void MaterialCohesive::computeNormal(const Array<Real> & position,
				     Array<Real> & normal,
				     GhostType ghost_type) {
}

/* -------------------------------------------------------------------------- */
inline UInt MaterialCohesive::getNbDataForElements(const Array<Element> & elements,
						   SynchronizationTag tag) const {

  switch (tag) {
  case _gst_smm_stress: {
    return 2 * spatial_dimension * sizeof(Real) * this->getModel().getNbQuadraturePoints(elements, "CohesiveFEEngine");
  }
  case _gst_smmc_damage: {
    return sizeof(Real) * this->getModel().getNbQuadraturePoints(elements, "CohesiveFEEngine");
  }
  default: {}
  }

  return 0;
}

/* -------------------------------------------------------------------------- */
inline void MaterialCohesive::packElementData(CommunicationBuffer & buffer,
					      const Array<Element> & elements,
					      SynchronizationTag tag) const {
  switch (tag) {
  case _gst_smm_stress: {
    packElementDataHelper(tractions, buffer, elements, "CohesiveFEEngine");
    packElementDataHelper(contact_tractions, buffer, elements, "CohesiveFEEngine");
    break;
  }
  case _gst_smmc_damage:
    packElementDataHelper(damage, buffer, elements, "CohesiveFEEngine"); break;
  default: {}
  }
}

/* -------------------------------------------------------------------------- */
inline void MaterialCohesive::unpackElementData(CommunicationBuffer & buffer,
						const Array<Element> & elements,
						SynchronizationTag tag) {
  switch (tag) {
  case _gst_smm_stress: {
    unpackElementDataHelper(tractions, buffer, elements, "CohesiveFEEngine");
    unpackElementDataHelper(contact_tractions, buffer, elements, "CohesiveFEEngine");
    break;
  }
  case _gst_smmc_damage:
    unpackElementDataHelper(damage, buffer, elements, "CohesiveFEEngine"); break;
  default: {}
  }
}
