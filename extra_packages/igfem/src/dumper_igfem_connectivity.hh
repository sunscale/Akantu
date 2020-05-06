/**
 * @file   dumper_igfem_connectivity.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Iterator for the IGFEM connectivity
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#ifndef __AKANTU_DUMPER_IGFEM_CONNECTIVITY_HH__
#define __AKANTU_DUMPER_IGFEM_CONNECTIVITY_HH__
/* -------------------------------------------------------------------------- */
#include "dumper_igfem_element_iterator.hh"
#include "dumper_igfem_generic_elemental_field.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
namespace dumpers {
/* -------------------------------------------------------------------------- */

template <class types>
class igfem_connectivity_field_iterator
    : public igfem_element_iterator<types, igfem_connectivity_field_iterator> {

public:
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

  typedef igfem_element_iterator<types,
                                 dumpers::igfem_connectivity_field_iterator>
      parent;
  typedef typename types::return_type return_type;
  typedef typename types::field_type field_type;
  typedef typename types::array_iterator array_iterator;

public:
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  igfem_connectivity_field_iterator(
      const field_type & field, const typename field_type::type_iterator & t_it,
      const typename field_type::type_iterator & t_it_end,
      const array_iterator & array_it, const array_iterator & array_it_end,
      const GhostType ghost_type = _not_ghost, UInt sub_element = 0)
      : parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type,
               sub_element) {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  return_type operator*() {
    const Vector<UInt> & element_connect = *this->array_it;

    /// get the local sub_element connectivity and the nodes per sub-element
    UInt * sub_connec_ptr =
        IGFEMHelper::getSubElementConnectivity(*this->tit, this->sub_element);
    UInt nb_nodes_sub_el =
        IGFEMHelper::getNbNodesPerSubElement(*this->tit, this->sub_element);

    /// get the global sub element connectivity
    Vector<UInt> sub_element_connect(nb_nodes_sub_el);
    for (UInt i = 0; i < nb_nodes_sub_el; ++i) {
      UInt lc = sub_connec_ptr[i];
      sub_element_connect(i) = element_connect(lc);
    }

    return sub_element_connect;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

private:
};

/* -------------------------------------------------------------------------- */
class IGFEMConnectivityField
    : public IGFEMGenericElementalField<SingleType<UInt, Vector, false>,
                                        igfem_connectivity_field_iterator> {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

public:
  typedef SingleType<UInt, Vector, false> types;
  typedef igfem_connectivity_field_iterator<types> iterator;
  typedef types::field_type field_type;
  typedef IGFEMGenericElementalField<types, igfem_connectivity_field_iterator>
      parent;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  IGFEMConnectivityField(const field_type & field,
                         UInt spatial_dimension = _all_dimensions,
                         GhostType ghost_type = _not_ghost)
      : parent(field, spatial_dimension, ghost_type) {}
};

/* -------------------------------------------------------------------------- */

} // namespace dumpers
} // namespace akantu

/* -------------------------------------------------------------------------- */
#endif /*__AKANTU_DUMPER_IGFEM_CONNECTIVITY_HH__ */
