/**
 * @file   petsc_matrix_inline_impl.cc
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Jul 21 21:42:12 2014
 *
 * @brief  Implementation of PetscMatrix inline functions
 *
 * @section LICENSE
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

// /* -------------------------------------------------------------------------- */
// inline UInt PetscMatrix::addToProfile(UInt i, UInt j) {

//   // initialize new value
//   MatSetValuesLocal(this->mat,1,i,1,j,0,INSERT_VALUES);

//   nb_non_zero++;

//   return nb_non_zero - 1;
// }

// /* -------------------------------------------------------------------------- */
// inline void PetscMatrix::clearProfile() {
//   //retain matrix but reinitialize its content
//   MatZeroEntries(this->mat);
//   nb_non_zero = 0;
// }

// /* -------------------------------------------------------------------------- */
// inline void PetscMatrix::clear() {
//   memset(a.storage(), 0, nb_non_zero*sizeof(Real));
// }

// /* -------------------------------------------------------------------------- */
// inline void PetscMatrix::addToMatrix(UInt i, UInt j, Real value) {
//   MatSetValue(this->mat,1,i,1,j,value,ADD_VALUES);
// }

// /* -------------------------------------------------------------------------- */
// inline Real PetscMatrix::operator()(UInt i, UInt j) const {
//   KeyCOO jcn_irn = key(i, j);
//   coordinate_list_map::const_iterator irn_jcn_k_it = irn_jcn_k.find(jcn_irn);
//   if(irn_jcn_k_it == irn_jcn_k.end()) return 0;
//   return this->mat(i,j).
// }

// /* -------------------------------------------------------------------------- */
// inline Real & PetscMatrix::operator()(UInt i, UInt j) {
//   KeyCOO jcn_irn = key(i, j);
//   coordinate_list_map::iterator irn_jcn_k_it = irn_jcn_k.find(jcn_irn);
//   AKANTU_DEBUG_ASSERT(irn_jcn_k_it != irn_jcn_k.end(),
// 		      "Couple (i,j) = (" << i << "," << j << ") does not exist in the profile");

//   return a.storage()[irn_jcn_k_it->second];
// }

