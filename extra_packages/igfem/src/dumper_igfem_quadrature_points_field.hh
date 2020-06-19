/**
 * @file   dumper_igfem_quadrature_points_field.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  description of IGFEM quadrature points field
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#ifndef __AKANTU_DUMPER_IGFEM_QUADRATURE_POINTS_FIELD_HH__
#define __AKANTU_DUMPER_IGFEM_QUADRATURE_POINTS_FIELD_HH__

/* -------------------------------------------------------------------------- */
#include "dumper_igfem_elemental_field.hh"

namespace akantu {
namespace dumpers {

/* -------------------------------------------------------------------------- */
template <typename types>
class igfem_quadrature_point_iterator
    : public igfem_element_iterator<types, igfem_quadrature_point_iterator> {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef igfem_element_iterator<types, dumpers::igfem_quadrature_point_iterator>
      parent;
  typedef typename types::data_type data_type;
  typedef typename types::return_type return_type;
  typedef typename types::field_type field_type;
  typedef typename types::array_iterator array_iterator;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  igfem_quadrature_point_iterator(
      const field_type & field, const typename field_type::type_iterator & t_it,
      const typename field_type::type_iterator & t_it_end,
      const array_iterator & array_it, const array_iterator & array_it_end,
      const GhostType ghost_type = _not_ghost, UInt sub_element = 0)
      : parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type,
               sub_element) {}

  return_type operator*() {
    const Vector<data_type> & mat_internal_field = *this->array_it;
    /// get nb data per sub element
    UInt nb_sub_1_internal_points =
        IGFEMHelper::getNbQuadraturePoints(*this->tit, 0);
    UInt nb_sub_2_internal_points =
        IGFEMHelper::getNbQuadraturePoints(*this->tit, 1);

    UInt nb_data = this->getNbDataPerElem(*(this->tit)) /
                   (nb_sub_1_internal_points + nb_sub_2_internal_points);

    UInt nb_sub_components = 0;
    if (!(this->sub_element))
      nb_sub_components = nb_data * nb_sub_1_internal_points;
    else
      nb_sub_components = nb_data * nb_sub_2_internal_points;

    Vector<data_type> sub_mat_internal_field(nb_sub_components);

    if (!(this->sub_element)) {
      for (UInt i = 0; i < nb_sub_components; ++i)
        sub_mat_internal_field(i) = mat_internal_field(i);
    } else {
      for (UInt i = 0; i < nb_sub_components; ++i)
        sub_mat_internal_field(i) =
            mat_internal_field(nb_data * nb_sub_1_internal_points + i);
    }

    return sub_mat_internal_field;
  }
};

// /* --------------------------------------------------------------------------
// */
// /* Fields type description */
// /* --------------------------------------------------------------------------
// */
// template<class types, template <class> class iterator_type>
// class GenericQuadraturePointsField :
//   public GenericElementalField<types,iterator_type> {

// public:

//   /* ------------------------------------------------------------------------
//   */
//   /* Typedefs */
//   /* ------------------------------------------------------------------------
//   */

//   typedef iterator_type<types> iterator;
//   typedef typename types::field_type field_type;
//   typedef typename iterator::it_type T;
//   typedef GenericElementalField<types,iterator_type> parent;

//   /* ------------------------------------------------------------------------
//   */
//   /* Constructors/Destructors */
//   /* ------------------------------------------------------------------------
//   */

//   GenericQuadraturePointsField(const field_type & field,
// 			       UInt spatial_dimension = _all_dimensions,
// 			       GhostType ghost_type = _not_ghost,
// 			       ElementKind element_kind = _ek_not_defined) :
//     parent(field, spatial_dimension, ghost_type, element_kind) { }

//   /* ------------------------------------------------------------------------
//   */
//   /* Methods */
//   /* ------------------------------------------------------------------------
//   */

//   virtual iterator begin() {
//     iterator it = parent::begin();
//     return it;
//   }

//   virtual iterator end  () {
//     iterator it = parent::end();
//     return it;
//   }

// };

/* -------------------------------------------------------------------------- */

} // namespace dumpers
} // namespace akantu

#endif /* __AKANTU_DUMPER_IGFEM_QUADRATURE_POINTS_FIELD_HH__ */
