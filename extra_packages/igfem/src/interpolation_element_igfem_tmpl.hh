/**
 * @file   element_class_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Interpolation property description for IGFEM
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
namespace akantu {
/* -------------------------------------------------------------------------- */
/* IGFEM elements                                                             */
/* -------------------------------------------------------------------------- */
#if defined(AKANTU_IGFEM)
/* -------------------------------------------------------------------------- */
template <> struct InterpolationProperty<_itp_igfem_segment_3> {
  static constexpr InterpolationKind kind{_itk_igfem};
  static constexpr UInt nb_nodes_per_element{3};
  static constexpr UInt natural_space_dimension{1};
  static constexpr InterpolationType parent_interpolation_type{
      _itp_lagrange_segment_2};
  static constexpr InterpolationType sub_iterpolation_type_1{
      _itp_lagrange_segment_2};
  static constexpr InterpolationType sub_interpolation_type_2{
      _itp_lagrange_segment_2};
  static constexpr UInt nb_sub_elements{2};
  static constexpr std::array<UInt, nb_sub_elements> sub_element_nb_nodes{
      {2, 2}};
  static constexpr std::array<UInt, 4> sub_element_connectivity_vect{
      {// first type
       0, 2,
       // second type
       2, 1}};
  static constexpr UInt nb_enrichments{1};
  static constexpr std::array<UInt, nb_enrichments * nb_sub_elements>
      enrichment_vect{{// on first subelement
                       1,
                       // on second subelement
                       0}};
};

/* -------------------------------------------------------------------------- */
template <> struct InterpolationProperty<_itp_igfem_triangle_4> {
  static constexpr InterpolationKind kind{_itk_igfem};
  static constexpr UInt nb_nodes_per_element{4};
  static constexpr UInt natural_space_dimension{2};
  static constexpr InterpolationType parent_interpolation_type{
      _itp_lagrange_triangle_3};
  static constexpr InterpolationType sub_iterpolation_type_1{
      _itp_lagrange_triangle_3};
  static constexpr InterpolationType sub_interpolation_type_2{
      _itp_lagrange_triangle_3};
  static constexpr UInt nb_sub_elements{2};
  static constexpr std::array<UInt, nb_sub_elements> sub_element_nb_nodes{
      {3, 3}};
  static constexpr std::array<UInt, 6> sub_element_connectivity_vect{
      {// irst type
       0, 1, 3,
       // second type
       0, 3, 2}};
  static constexpr UInt nb_enrichments{1};
  static constexpr std::array<UInt, nb_enrichments * nb_sub_elements>
      enrichment_vect{{// on first subelement
                       2,
                       // on second subelement
                       1}};
};

/* -------------------------------------------------------------------------- */
template <> struct InterpolationProperty<_itp_igfem_triangle_5> {
  static constexpr InterpolationKind kind{_itk_igfem};
  static constexpr UInt nb_nodes_per_element{5};
  static constexpr UInt natural_space_dimension{2};
  static constexpr InterpolationType parent_interpolation_type{
      _itp_lagrange_triangle_3};
  static constexpr InterpolationType sub_iterpolation_type_1{
      _itp_lagrange_triangle_3};
  static constexpr InterpolationType sub_interpolation_type_2{
      _itp_lagrange_quadrangle_4};
  static constexpr UInt nb_sub_elements{2};
  static constexpr std::array<UInt, nb_sub_elements> sub_element_nb_nodes{
      {3, 4}};
  static constexpr std::array<UInt, 7> sub_element_connectivity_vect{
      {// first type
       0, 3, 4,
       // second type
       3, 1, 2, 4}};
  static constexpr UInt nb_enrichments{2};
  static constexpr std::array<UInt, nb_enrichments * nb_sub_elements>
      enrichment_vect{{// on first subelement
                       1, 2,
                       // on second subelement
                       0, 3}};
};
#endif
} // namespace akantu
