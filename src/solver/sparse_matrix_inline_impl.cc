/**
 * @file   sparse_matrix_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  implementation of inline methods of the SparseMatrix class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
inline UInt SparseMatrix::addToProfile(UInt i, UInt j) {

  KeyCOO jcn_irn = key(i, j);

  coordinate_list_map::iterator it = irn_jcn_k.find(jcn_irn);
  if (!(it == irn_jcn_k.end()))
      return it->second;

//  AKANTU_DEBUG_ASSERT(irn_jcn_k.find(jcn_irn) == irn_jcn_k.end(),
//		      "Couple (i,j) = (" << i << "," << j << ") already in the profile");

  irn.push_back(i + 1);
  jcn.push_back(j + 1);

  //  a.resize(a.getSize() + 1);
  Real zero = 0;
  a.push_back(zero);

  irn_jcn_k[jcn_irn] = nb_non_zero;

  nb_non_zero++;

  return nb_non_zero - 1;
}

/* -------------------------------------------------------------------------- */
inline void SparseMatrix::clearProfile() {
  irn_jcn_k.clear();
  irn.resize(0);
  jcn.resize(0);
  a.resize(0);
  nb_non_zero = 0;
}

/* -------------------------------------------------------------------------- */
inline void SparseMatrix::clear() {
  memset(a.storage(), 0, nb_non_zero*sizeof(Real));
}

/* -------------------------------------------------------------------------- */
inline void SparseMatrix::addToMatrix(UInt i, UInt j, Real value) {
  KeyCOO jcn_irn = key(i, j);
  coordinate_list_map::iterator irn_jcn_k_it = irn_jcn_k.find(jcn_irn);

  AKANTU_DEBUG_ASSERT(irn_jcn_k_it != irn_jcn_k.end(),
		      "Couple (i,j) = (" << i << "," << j << ") does not exist in the profile");

  a.storage()[irn_jcn_k_it->second] += value;
}

/* -------------------------------------------------------------------------- */
inline Real SparseMatrix::operator()(UInt i, UInt j) const {
  KeyCOO jcn_irn = key(i, j);
  coordinate_list_map::const_iterator irn_jcn_k_it = irn_jcn_k.find(jcn_irn);
  if(irn_jcn_k_it == irn_jcn_k.end()) return 0;
  return a.storage()[irn_jcn_k_it->second];
}

/* -------------------------------------------------------------------------- */
inline Real & SparseMatrix::operator()(UInt i, UInt j) {
  KeyCOO jcn_irn = key(i, j);
  coordinate_list_map::iterator irn_jcn_k_it = irn_jcn_k.find(jcn_irn);
  AKANTU_DEBUG_ASSERT(irn_jcn_k_it != irn_jcn_k.end(),
		      "Couple (i,j) = (" << i << "," << j << ") does not exist in the profile");

  return a.storage()[irn_jcn_k_it->second];
}

/* -------------------------------------------------------------------------- */
// inline void SparseMatrix::addToMatrixSym(const Array<Real> & local_matrix,
// 					 const Element & element) {
//   AKANTU_DEBUG_ASSERT(element_to_sparse_profile[element.type] != NULL,
// 		      "No profile stored for this kind of element call first buildProfile()");

//   UInt nb_values_per_elem   = element_to_sparse_profile[element.type]->getNbComponent();
//   UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

//   Real * mat_val = local_matrix.storage();
//   UInt * elem_to_sparse_val = element_to_sparse_profile[element.type]->values + element.element * nb_values_per_elem;
//   Real * a_val = a.storage();

//   for (UInt j = 0; j < nb_nodes_per_element * nb_degree_of_freedom; ++j) {
//     UInt i_end = (sparse_matrix_type == _symmetric) ? j + 1 :
//       nb_nodes_per_element * nb_degree_of_freedom;
//     for (UInt i = 0; i < i_end; ++i) {
//       UInt k = *(elem_to_sparse_val++);
//       a_val[k] += *(mat_val++);
//     }
//   }
// }

/* -------------------------------------------------------------------------- */
// inline void SparseMatrix::addToMatrix(Real * local_matrix,
// 				      const Element & element,
// 				      UInt nb_nodes_per_element) {
//   AKANTU_DEBUG_ASSERT(element_to_sparse_profile[element.type] != NULL,
// 		      "No profile stored for this kind of element call first buildProfile()");

//   UInt nb_values_per_elem   = element_to_sparse_profile[element.type]->getNbComponent();
//   //  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

//   Real * mat_val = local_matrix;
//   UInt * elem_to_sparse_val = element_to_sparse_profile[element.type]->values + element.element * nb_values_per_elem;
//   Real * a_val = a.storage();

//   for (UInt i = 0; i < nb_nodes_per_element * nb_degree_of_freedom; ++i) {
//     UInt j_start = (sparse_matrix_type == _symmetric) ? i : 0;
//     UInt elem_to_sparse_i = i * nb_nodes_per_element * nb_degree_of_freedom;
//     for (UInt j = j_start; j < nb_nodes_per_element * nb_degree_of_freedom; ++j) {
//       UInt k = elem_to_sparse_val[elem_to_sparse_i + j];
//       a_val[k] += mat_val[elem_to_sparse_i + j];
//     }
//   }
// }
