/**
 * @file   material_igfem_saw_tooth_damage.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Implementation of the squentially linear saw-tooth damage model for
 * IGFEM elements
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */
/* -------------------------------------------------------------------------- */
#include "material_igfem_saw_tooth_damage.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
MaterialIGFEMSawToothDamage<dim>::MaterialIGFEMSawToothDamage(
    SolidMechanicsModel & model, const ID & id)
    : Material(model, id), Parent(model, id), Sc("Sc", *this),
      equivalent_stress("equivalent_stress", *this),
      norm_max_equivalent_stress(0) {
  AKANTU_DEBUG_IN();

  this->Sc.initialize(1);
  this->equivalent_stress.initialize(1);
  this->damage.setElementKind(_ek_igfem);
  this->damage.setFEEngine(*(this->fem));

  this->internals_to_transfer.push_back("damage");
  this->internals_to_transfer.push_back("Sc");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim> void MaterialIGFEMSawToothDamage<dim>::initMaterial() {
  AKANTU_DEBUG_IN();

  Parent::initMaterial();

  /// get the parameters for the sub-material that can be damaged
  ID mat_name = this->sub_material_names[1];
  const Material & mat = this->model->getMaterial(mat_name);
  this->prescribed_dam = mat.getParam<Real>("prescribed_dam");
  this->dam_threshold = mat.getParam<Real>("dam_threshold");
  this->dam_tolerance = mat.getParam<Real>("dam_tolerance");
  this->max_damage = mat.getParam<Real>("max_damage");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialIGFEMSawToothDamage<spatial_dimension>::
    computeNormalizedEquivalentStress(const Array<Real> & grad_u,
                                      ElementType el_type,
                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// Vector to store eigenvalues of current stress tensor
  Vector<Real> eigenvalues(spatial_dimension);

  Array<Real>::const_iterator<Real> Sc_it = Sc(el_type, ghost_type).begin();
  Array<Real>::iterator<Real> equivalent_stress_it =
      equivalent_stress(el_type, ghost_type).begin();

  Array<Real>::const_matrix_iterator grad_u_it =
      grad_u.begin(spatial_dimension, spatial_dimension);
  Array<Real>::const_matrix_iterator grad_u_end =
      grad_u.end(spatial_dimension, spatial_dimension);
  /// get pointer to internals
  Real * lambda_ptr = this->lambda(el_type, ghost_type).storage();
  Real * mu_ptr = this->mu(el_type, ghost_type).storage();
  Real * dam = this->damage(el_type, ghost_type).storage();
  Matrix<Real> sigma(spatial_dimension, spatial_dimension);
  for (; grad_u_it != grad_u_end; ++grad_u_it) {
    sigma.clear();
    MaterialIGFEMElastic<spatial_dimension>::computeStressOnQuad(
        *grad_u_it, sigma, *lambda_ptr, *mu_ptr);
    computeDamageAndStressOnQuad(sigma, *dam);

    /// compute eigenvalues
    sigma.eig(eigenvalues);
    /// find max eigenvalue and normalize by tensile strength
    *equivalent_stress_it =
        *(std::max_element(eigenvalues.storage(),
                           eigenvalues.storage() + spatial_dimension)) /
        *(Sc_it);
    ++Sc_it;
    ++equivalent_stress_it;
    ++dam;
    ++lambda_ptr;
    ++mu_ptr;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialIGFEMSawToothDamage<spatial_dimension>::computeAllStresses(
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// reset normalized maximum equivalent stress
  if (ghost_type == _not_ghost)
    norm_max_equivalent_stress = 0;

  Parent::computeAllStresses(ghost_type);

  /// find global Gauss point with highest stress
  StaticCommunicator & comm =
      akantu::StaticCommunicator::getStaticCommunicator();
  comm.allReduce(&norm_max_equivalent_stress, 1, _so_max);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialIGFEMSawToothDamage<spatial_dimension>::
    findMaxNormalizedEquivalentStress(ElementType el_type,
                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if (ghost_type == _not_ghost) {

    const Array<Real> & e_stress = equivalent_stress(el_type);
    Array<Real>::const_iterator<Real> equivalent_stress_it = e_stress.begin();
    Array<Real>::const_iterator<Real> equivalent_stress_end = e_stress.end();
    Array<Real> & dam = this->damage(el_type);
    Array<Real>::iterator<Real> dam_it = dam.begin();
    Array<UInt> & sub_mat = this->sub_material(el_type);
    Array<UInt>::iterator<UInt> sub_mat_it = sub_mat.begin();

    for (; equivalent_stress_it != equivalent_stress_end;
         ++equivalent_stress_it, ++dam_it, ++sub_mat_it) {
      /// check if max equivalent stress for this element type is greater than
      /// the current norm_max_eq_stress and if the element is not already fully
      /// damaged
      if ((*equivalent_stress_it > norm_max_equivalent_stress &&
           *dam_it < max_damage) &&
          *sub_mat_it != 0) {
        norm_max_equivalent_stress = *equivalent_stress_it;
      }
    }
  }
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialIGFEMSawToothDamage<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Parent::computeStress(el_type, ghost_type);

  Real * dam = this->damage(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computeDamageAndStressOnQuad(sigma, *dam);
  ++dam;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  computeNormalizedEquivalentStress(this->gradu(el_type, ghost_type), el_type,
                                    ghost_type);
  norm_max_equivalent_stress = 0;
  findMaxNormalizedEquivalentStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
UInt MaterialIGFEMSawToothDamage<spatial_dimension>::updateDamage() {
  UInt nb_damaged_elements = 0;
  AKANTU_DEBUG_ASSERT(prescribed_dam > 0.,
                      "Your prescribed damage must be greater than zero");

  if (norm_max_equivalent_stress >= 1.) {

    AKANTU_DEBUG_IN();
    GhostType ghost_type = _not_ghost;
    ;

    Mesh::type_iterator it = this->model->getMesh().firstType(
        spatial_dimension, ghost_type, _ek_igfem);
    Mesh::type_iterator last_type = this->model->getMesh().lastType(
        spatial_dimension, ghost_type, _ek_igfem);

    for (; it != last_type; ++it) {
      ElementType el_type = *it;
      Array<UInt> & sub_mat = this->sub_material(el_type);
      Array<UInt>::iterator<UInt> sub_mat_it = sub_mat.begin();
      const Array<Real> & e_stress = equivalent_stress(el_type);
      Array<Real>::const_iterator<Real> equivalent_stress_it = e_stress.begin();
      Array<Real>::const_iterator<Real> equivalent_stress_end = e_stress.end();
      Array<Real> & dam = this->damage(el_type);
      Array<Real>::iterator<Real> dam_it = dam.begin();

      /// loop over all the elements of the given type
      UInt nb_element = this->element_filter(el_type, ghost_type).getSize();
      UInt nb_quads = this->fem->getNbIntegrationPoints(el_type, ghost_type);
      bool damage_element = false;
      for (UInt e = 0; e < nb_element; ++e) {
        damage_element = false;
        /// check if damage occurs in the element
        for (UInt q = 0; q < nb_quads; ++q) {
          if (*equivalent_stress_it >=
                  (1 - dam_tolerance) * norm_max_equivalent_stress &&
              *sub_mat_it != 0) {
            damage_element = true;
          }
          ++sub_mat_it;
          ++equivalent_stress_it;
        }

        if (damage_element) {
          nb_damaged_elements += 1;
          /// damage the element
          sub_mat_it -= nb_quads;
          for (UInt q = 0; q < nb_quads; ++q) {
            if (*sub_mat_it) {
              if (*dam_it < dam_threshold)
                *dam_it += prescribed_dam;
              else
                *dam_it = max_damage;
            }
            ++sub_mat_it;
            ++dam_it;
          }
        } else {
          dam_it += nb_quads;
        }
      }
    }
  }

  StaticCommunicator & comm =
      akantu::StaticCommunicator::getStaticCommunicator();
  comm.allReduce(&nb_damaged_elements, 1, _so_sum);
  AKANTU_DEBUG_OUT();
  return nb_damaged_elements;
}
/* -------------------------------------------------------------------------- */
// template<UInt spatial_dimension>
// void
// MaterialIGFEMSawToothDamage<spatial_dimension>::updateEnergiesAfterDamage(ElementType
// el_type, GhostType ghost_type) {
//   MaterialDamage<spatial_dimension>::updateEnergies(el_type, ghost_type);
// }

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialIGFEMSawToothDamage<spatial_dimension>::onElementsAdded(
    __attribute__((unused)) const Array<Element> & element_list,
    __attribute__((unused)) const NewElementsEvent & event) {

  /// update elastic constants
  MaterialIGFEMElastic<spatial_dimension>::onElementsAdded(element_list, event);
  IGFEMInternalField<Real> quadrature_points_coordinates(
      "quadrature_points_coordinates", *this);
  quadrature_points_coordinates.initialize(spatial_dimension);
  this->computeQuadraturePointsCoordinates(quadrature_points_coordinates,
                                           _not_ghost);
  this->computeQuadraturePointsCoordinates(quadrature_points_coordinates,
                                           _ghost);

  Array<Element>::const_iterator<Element> el_begin = element_list.begin();
  Array<Element>::const_iterator<Element> el_end = element_list.end();
  Element el;
  el.kind = _ek_igfem;
  /// loop over all the elements in the filter
  for (ghost_type_t::iterator g = ghost_type_t::begin();
       g != ghost_type_t::end(); ++g) {
    GhostType ghost_type = *g;
    /// loop over all types in the material
    typedef ElementTypeMapArray<UInt>::type_iterator iterator;
    iterator it = this->element_filter.firstType(spatial_dimension, ghost_type,
                                                 _ek_igfem);
    iterator last_type =
        this->element_filter.lastType(spatial_dimension, ghost_type, _ek_igfem);
    for (; it != last_type; ++it) {
      /// store the elements added to this material
      std::vector<Element> new_elements;
      Array<UInt> & elem_filter = this->element_filter(*it, ghost_type);
      UInt nb_element = elem_filter.getSize();
      UInt nb_quads =
          quadrature_points_coordinates.getFEEngine().getNbIntegrationPoints(
              *it, ghost_type);
      Array<Real> added_quads(0, spatial_dimension);
      const Array<Real> & quads =
          quadrature_points_coordinates(*it, ghost_type);
      Array<Real>::const_vector_iterator quad = quads.begin(spatial_dimension);
      el.type = *it;
      el.ghost_type = ghost_type;

      for (UInt e = 0; e < nb_element; ++e) {
        /// global number of element
        el.element = elem_filter(e);
        if (std::find(el_begin, el_end, el) != el_end) {
          new_elements.push_back(el);
          for (UInt q = 0; q < nb_quads; ++q) {
            added_quads.push_back(*quad);
            ++quad;
          }
        } else
          quad += nb_quads;
      }

      /// get the extrapolated values
      for (UInt i = 0; i < this->internals_to_transfer.size(); ++i) {
        if (!new_elements.size())
          continue;
        const ID name = this->getID() + ":" + this->internals_to_transfer[i];
        Array<Real> & internal =
            (*(this->internal_vectors_real[name]))(*it, ghost_type);
        UInt nb_component = internal.getNbComponent();
        /// Array<Real>::vector_iterator internal_it =
        /// internal.begin(nb_component);
        Array<Real> extrapolated_internal_values(new_elements.size() * nb_quads,
                                                 nb_component);
        this->model->transferInternalValues(this->internals_to_transfer[i],
                                            new_elements, added_quads,
                                            extrapolated_internal_values);

        UInt * sub_mat_ptr = this->sub_material(*it, ghost_type).storage();
        Array<Real>::const_vector_iterator extrapolated_it =
            extrapolated_internal_values.begin(nb_component);
        Array<Real>::vector_iterator internal_it = internal.begin(nb_component);
        /// copy back the results
        for (UInt j = 0; j < new_elements.size(); ++j) {
          Element element_local = this->convertToLocalElement(new_elements[j]);
          for (UInt q = 0; q < nb_quads; ++q, ++extrapolated_it) {
            if (sub_mat_ptr[element_local.element * nb_quads + q])
              internal_it[element_local.element * nb_quads + q] =
                  *extrapolated_it;
          }
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
// template<UInt spatial_dimension>
// void
// MaterialIGFEMSawToothDamage<spatial_dimension>::transferInternals(Material &
// old_mat,
// 								       std::vector<ElementPair> & element_pairs)
// {

// Element new_el_global;
// Element new_el_local;
// Element old_el_global;
// Element old_el_local;

// /// get the fe-engine of the old material
// FEEngine & fem_old_mat = old_mat.getFEEngine();
// for (UInt e = 0; e < element_pairs.size(); ++e) {
//   new_el_global = element_pairs[e].first;
//   old_el_global = element_pairs[e].second;
//   /// get the number of the elements in their materials
//   old_el_local = old_el_global;
//   Array<UInt> & mat_local_numbering =
//   this->model->getMaterialLocalNumbering(old_el_global.type,
//   old_el_global.ghost_type);
//   old_el_local.element = mat_local_numbering(old_el_global.element);
//   new_el_local = this->convertToLocalElement(new_el_global);
//   UInt nb_old_quads = fem_old_mat.getNbIntegrationPoints(old_el_global.type,
//   old_el_global.ghost_type);
//   UInt nb_new_quads = this->fem->getNbIntegrationPoints(new_el_global.type,
//   new_el_global.ghost_type);

//   if (old_mat.isInternal("damage", Mesh::getKind(old_el_global.type))
// 	&& old_mat.isInternal("Sc", Mesh::getKind(old_el_global.type)) ) {

//     UInt quad;
//     Vector<Real> el_old_damage(nb_old_quads);
//     Vector<Real> el_old_Sc(nb_old_quads);
//     Vector<Real> el_new_damage(nb_new_quads);
//     Vector<Real> el_new_Sc(nb_new_quads);
//     const Array<Real> & old_Sc = old_mat.getArray("Sc", old_el_global.type,
//     old_el_global.ghost_type);
//     const Array<Real> & old_damage = old_mat.getArray("damage",
//     old_el_global.type, old_el_global.ghost_type);
//     Array<Real> & new_Sc = this->Sc(new_el_global.type,
//     new_el_global.ghost_type);
//     Array<Real> & new_damage = this->damage(new_el_global.type,
//     new_el_global.ghost_type);

//     for (UInt q = 0; q < nb_old_quads; ++q) {
// 	quad = old_el_local.element * nb_old_quads + q;
// 	el_old_damage(q) = old_damage(quad);
// 	el_old_Sc(q) = old_Sc(quad);
//     }

//     this->interpolateInternal(new_el_global, old_el_global, el_new_damage,
//     el_old_damage, nb_new_quads, nb_old_quads);
//     this->interpolateInternal(new_el_global, old_el_global, el_new_Sc,
//     el_old_Sc, nb_new_quads, nb_old_quads);

//     for (UInt q = 0; q < nb_new_quads; ++q) {
// 	quad = new_el_local.element * nb_new_quads + q;
// 	if (this->sub_material(new_el_global.type,new_el_global.ghost_type)(quad)) {
// 	  new_damage(quad) = el_new_damage(q);
// 	  new_Sc(quad) = el_new_damage(q);
// 	}
//     }
//   }

//   else
//     AKANTU_DEBUG_ASSERT((!old_mat.isInternal("damage",
//     Mesh::getKind(old_el_global.type))
// 			   && !old_mat.isInternal("Sc", Mesh::getKind(old_el_global.type))
// ),
// 			  "old material has damage or Sc but not both!!!!");
// }
// }

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(MaterialIGFEMSawToothDamage);

} // namespace akantu

/*     /// get the material index in the model from this current material
    UInt this_mat_idx =  this->model->getMaterialIndex(this->name);

    /// number of elements in this event
    UInt nb_new_elements = element_list.getSize();

    // const Array<Element> & old_elements = igfem_event->getOldElementsList();

    // std::map<UInt, std::vector<ElementPair> > elements_by_old_mat;

    /// store the old elements sorted by their material
    for (UInt e = 0; e < nb_new_elements; ++e) {
      const Element new_el = element_list(e);
      const Array<UInt> & mat_idx =
   this->model->getMaterialByElement(new_el.type, new_el.ghost_type);
      if ( mat_idx(new_el.element) != this_mat_idx )
    continue; /// new element is not part of this material: nothing to be done
      /// get the corresponding old element and store new and old one as pair
      const Element old_el = old_elements(e);
      ElementPair el_pair(new_el, old_el);
      const Array<UInt> & old_mat_idx =
   this->model->getMaterialByElement(old_el.type, old_el.ghost_type);
      UInt mat_old_idx = old_mat_idx(old_el.element);
      elements_by_old_mat[mat_old_idx].push_back(el_pair);
    }


    /// loop over all the element pairs in the map
    for (std::map<UInt, std::vector<ElementPair> >::iterator map_it =
   elements_by_old_mat.begin();
     map_it != elements_by_old_mat.end(); ++map_it) {
      /// get the vector of old and new element pairs
      std::vector<ElementPair > & element_pairs = map_it->second;
      Material & old_mat = this->model->getMaterial(map_it->first);
      this->transferInternals(old_mat, element_pairs);
    }
*/
