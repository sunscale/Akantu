/**
 * @file   material_cohesive_bilinear.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Feb 22 2012
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Bilinear cohesive constitutive law
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
#include "material_cohesive_bilinear.hh"
//#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialCohesiveBilinear<spatial_dimension>::MaterialCohesiveBilinear(
    SolidMechanicsModel & model, const ID & id)
    : MaterialCohesiveLinear<spatial_dimension>(model, id) {
  AKANTU_DEBUG_IN();

  this->registerParam("delta_0", delta_0, Real(0.),
                      _pat_parsable | _pat_readable,
                      "Elastic limit displacement");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  this->sigma_c_eff.setRandomDistribution(this->sigma_c.getRandomParameter());
  MaterialCohesiveLinear<spatial_dimension>::initMaterial();

  this->delta_max.setDefaultValue(delta_0);
  this->insertion_stress.setDefaultValue(0);

  this->delta_max.reset();
  this->insertion_stress.reset();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::onElementsAdded(
    const Array<Element> & element_list, const NewElementsEvent & event) {
  AKANTU_DEBUG_IN();

  MaterialCohesiveLinear<spatial_dimension>::onElementsAdded(element_list,
                                                             event);

  bool scale_traction = false;

  // don't scale sigma_c if volume_s hasn't been specified by the user
  if (!Math::are_float_equal(this->volume_s, 0.)) {
    scale_traction = true;
  }

  Array<Element>::const_scalar_iterator el_it = element_list.begin();
  Array<Element>::const_scalar_iterator el_end = element_list.end();

  for (; el_it != el_end; ++el_it) {
    // filter not ghost cohesive elements
    if ((el_it->ghost_type != _not_ghost) or
        (Mesh::getKind(el_it->type) != _ek_cohesive)) {
      continue;
    }

    UInt index = el_it->element;
    ElementType type = el_it->type;
    UInt nb_element = this->model->getMesh().getNbElement(type);
    UInt nb_quad_per_element = this->fem_cohesive.getNbIntegrationPoints(type);

    auto sigma_c_begin = this->sigma_c_eff(type).begin_reinterpret(
        nb_quad_per_element, nb_element);
    Vector<Real> sigma_c_vec = sigma_c_begin[index];

    auto delta_c_begin = this->delta_c_eff(type).begin_reinterpret(
        nb_quad_per_element, nb_element);
    Vector<Real> delta_c_vec = delta_c_begin[index];

    if (scale_traction) {
      scaleTraction(*el_it, sigma_c_vec);
    }

    /**
     * Recompute sigma_c as
     * @f$ {\sigma_c}_\textup{new} =
     * \frac{{\sigma_c}_\textup{old} \delta_c} {\delta_c - \delta_0} @f$
     */

    for (UInt q = 0; q < nb_quad_per_element; ++q) {
      delta_c_vec(q) = 2 * this->G_c / sigma_c_vec(q);

      if (delta_c_vec(q) - delta_0 < Math::getTolerance()) {
        AKANTU_ERROR("delta_0 = " << delta_0 << " must be lower than delta_c = "
                                  << delta_c_vec(q)
                                  << ", modify your material file");
      }

      sigma_c_vec(q) *= delta_c_vec(q) / (delta_c_vec(q) - delta_0);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::scaleTraction(
    const Element & el, Vector<Real> & sigma_c_vec) {
  AKANTU_DEBUG_IN();

  Real base_sigma_c = this->sigma_c_eff;

  const Mesh & mesh_facets = this->model->getMeshFacets();
  const FEEngine & fe_engine = this->model->getFEEngine();

  auto coh_element_to_facet_begin =
      mesh_facets.getSubelementToElement(el.type).begin(2);
  const Vector<Element> & coh_element_to_facet =
      coh_element_to_facet_begin[el.element];

  // compute bounding volume
  Real volume = 0;

  // loop over facets
  for (UInt f = 0; f < 2; ++f) {
    const Element & facet = coh_element_to_facet(f);

    const Array<std::vector<Element>> & facet_to_element =
        mesh_facets.getElementToSubelement(facet.type, facet.ghost_type);

    const std::vector<Element> & element_list = facet_to_element(facet.element);

    auto elem = element_list.begin();
    auto elem_end = element_list.end();

    // loop over elements connected to each facet
    for (; elem != elem_end; ++elem) {
      // skip cohesive elements and dummy elements
      if (*elem == ElementNull || Mesh::getKind(elem->type) == _ek_cohesive) {
        continue;
      }

      // unit vector for integration in order to obtain the volume
      UInt nb_quadrature_points = fe_engine.getNbIntegrationPoints(elem->type);
      Vector<Real> unit_vector(nb_quadrature_points, 1);

      volume += fe_engine.integrate(unit_vector, elem->type, elem->element,
                                    elem->ghost_type);
    }
  }

  // scale sigma_c
  sigma_c_vec -= base_sigma_c;
  sigma_c_vec *= std::pow(this->volume_s / volume, 1. / this->m_s);
  sigma_c_vec += base_sigma_c;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::computeTraction(
    const Array<Real> & normal, ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  MaterialCohesiveLinear<spatial_dimension>::computeTraction(normal, el_type,
                                                             ghost_type);

  // adjust damage
  auto delta_c_it = this->delta_c_eff(el_type, ghost_type).begin();
  auto delta_max_it = this->delta_max(el_type, ghost_type).begin();
  auto damage_it = this->damage(el_type, ghost_type).begin();
  auto damage_end = this->damage(el_type, ghost_type).end();

  for (; damage_it != damage_end; ++damage_it, ++delta_max_it, ++delta_c_it) {
    *damage_it =
        std::max((*delta_max_it - delta_0) / (*delta_c_it - delta_0), Real(0.));
    *damage_it = std::min(*damage_it, Real(1.));
  }
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(cohesive_bilinear, MaterialCohesiveBilinear);

} // namespace akantu
