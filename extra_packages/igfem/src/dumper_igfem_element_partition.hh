/**
 * @file   dumper_igfem_element_partition.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Element partition field for IGFEM sub-elements
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
namespace akantu {
namespace dumpers {

/* -------------------------------------------------------------------------- */
template <class types>
class igfem_element_partition_field_iterator
    : public igfem_element_iterator<types,
                                    igfem_element_partition_field_iterator> {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef igfem_element_iterator<types,
                                 dumpers::igfem_element_partition_field_iterator>
      parent;
  typedef typename types::return_type return_type;
  typedef typename types::array_iterator array_iterator;
  typedef typename types::field_type field_type;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  igfem_element_partition_field_iterator(
      const field_type & field, const typename field_type::type_iterator & t_it,
      const typename field_type::type_iterator & t_it_end,
      const array_iterator & array_it, const array_iterator & array_it_end,
      const GhostType ghost_type = _not_ghost, UInt sub_element = 0)
      : parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type,
               sub_element) {
    prank = StaticCommunicator::getStaticCommunicator().whoAmI();
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  return_type operator*() { return return_type(1, prank); }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  UInt prank;
};

/* -------------------------------------------------------------------------- */
template <bool filtered = false>
class IGFEMElementPartitionField : public IGFEMGenericElementalField<
                                       SingleType<UInt, Vector, filtered>,
                                       igfem_element_partition_field_iterator> {
public:
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

  typedef SingleType<UInt, Vector, filtered> types;
  typedef igfem_element_partition_field_iterator<types> iterator;
  typedef IGFEMGenericElementalField<types,
                                     igfem_element_partition_field_iterator>
      parent;
  typedef typename types::field_type field_type;

public:
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  IGFEMElementPartitionField(const field_type & field,
                             UInt spatial_dimension = _all_dimensions,
                             GhostType ghost_type = _not_ghost,
                             ElementKind kind = _ek_igfem)
      : parent(field, spatial_dimension, ghost_type, kind) {
    this->homogeneous = true;
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  UInt getDim() { return 1; }
};

/* -------------------------------------------------------------------------- */

} // namespace dumpers
} // namespace akantu
