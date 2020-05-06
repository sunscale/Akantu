/**
 * @file   igfem_helper.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Fri Jul  3 15:49:55 2015
 *
 * @brief
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "igfem_helper.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

template <>
ElementType ElementTypeIGFEMData<_not_defined>::igfem_element_types[] = {
    _not_defined};
template <> UInt ElementTypeIGFEMData<_not_defined>::nb_igfem_types = 0;
/* -------------------------------------------------------------------------- */
template <>
ElementType ElementTypeIGFEMData<_point_1>::igfem_element_types[] = {
    _not_defined};
template <> UInt ElementTypeIGFEMData<_point_1>::nb_igfem_types = 0;
/* -------------------------------------------------------------------------- */
template <>
ElementType ElementTypeIGFEMData<_segment_2>::igfem_element_types[] = {
    _igfem_segment_3};
template <> UInt ElementTypeIGFEMData<_segment_2>::nb_igfem_types = 1;
/* -------------------------------------------------------------------------- */
template <>
ElementType ElementTypeIGFEMData<_segment_3>::igfem_element_types[] = {
    _not_defined};
template <> UInt ElementTypeIGFEMData<_segment_3>::nb_igfem_types = 0;
/* -------------------------------------------------------------------------- */
template <>
ElementType ElementTypeIGFEMData<_triangle_3>::igfem_element_types[] = {
    _igfem_triangle_4, _igfem_triangle_5};
template <> UInt ElementTypeIGFEMData<_triangle_3>::nb_igfem_types = 2;
/* -------------------------------------------------------------------------- */
template <>
ElementType ElementTypeIGFEMData<_triangle_6>::igfem_element_types[] = {
    _not_defined};
template <> UInt ElementTypeIGFEMData<_triangle_6>::nb_igfem_types = 0;
/* -------------------------------------------------------------------------- */
template <>
ElementType ElementTypeIGFEMData<_tetrahedron_4>::igfem_element_types[] = {
    _not_defined};
template <> UInt ElementTypeIGFEMData<_tetrahedron_4>::nb_igfem_types = 0;
/* -------------------------------------------------------------------------- */
template <>
ElementType ElementTypeIGFEMData<_tetrahedron_10>::igfem_element_types[] = {
    _not_defined};
template <> UInt ElementTypeIGFEMData<_tetrahedron_10>::nb_igfem_types = 0;
/* -------------------------------------------------------------------------- */
template <>
ElementType ElementTypeIGFEMData<_quadrangle_4>::igfem_element_types[] = {
    _not_defined};
template <> UInt ElementTypeIGFEMData<_quadrangle_4>::nb_igfem_types = 0;
/* -------------------------------------------------------------------------- */
template <>
ElementType ElementTypeIGFEMData<_quadrangle_8>::igfem_element_types[] = {
    _not_defined};
template <> UInt ElementTypeIGFEMData<_quadrangle_8>::nb_igfem_types = 0;
/* -------------------------------------------------------------------------- */
template <>
ElementType ElementTypeIGFEMData<_hexahedron_8>::igfem_element_types[] = {
    _not_defined};
template <> UInt ElementTypeIGFEMData<_hexahedron_8>::nb_igfem_types = 0;
/* -------------------------------------------------------------------------- */
template <>
ElementType ElementTypeIGFEMData<_hexahedron_20>::igfem_element_types[] = {
    _not_defined};
template <> UInt ElementTypeIGFEMData<_hexahedron_20>::nb_igfem_types = 0;
/* -------------------------------------------------------------------------- */
template <>
ElementType ElementTypeIGFEMData<_pentahedron_6>::igfem_element_types[] = {
    _not_defined};
template <> UInt ElementTypeIGFEMData<_pentahedron_6>::nb_igfem_types = 0;
/* -------------------------------------------------------------------------- */
template <>
ElementType ElementTypeIGFEMData<_pentahedron_15>::igfem_element_types[] = {
    _not_defined};
template <> UInt ElementTypeIGFEMData<_pentahedron_15>::nb_igfem_types = 0;
/* -------------------------------------------------------------------------- */

} // namespace akantu
