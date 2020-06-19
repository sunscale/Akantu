/**
 * @file   sparse_solver_inline_impl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Sun Aug 13 2017
 *
 * @brief  implementation of solver inline functions
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
// inline UInt Solver::getNbDataForDOFs(const Array<UInt> & dofs,
// 				     SynchronizationTag tag) const {
//   AKANTU_DEBUG_IN();

//   UInt size = 0;

//   switch(tag) {
//   case SynchronizationTag::_solver_solution: {
//     size += dofs.size() * sizeof(Real);
//     break;
//   }
//   default: {  }
//   }

//   AKANTU_DEBUG_OUT();
//   return size;
// }

// /* --------------------------------------------------------------------------
// */
// inline void Solver::packDOFData(CommunicationBuffer & buffer,
// 				const Array<UInt> & dofs,
// 				SynchronizationTag tag) const {
//   AKANTU_DEBUG_IN();

//   switch(tag) {
//   case SynchronizationTag::_solver_solution: {
//     packDOFDataHelper(*solution, buffer, dofs);
//     break;
//   }
//   default: {
//   }
//   }

//   AKANTU_DEBUG_OUT();
// }

// /* --------------------------------------------------------------------------
// */
// inline void Solver::unpackDOFData(CommunicationBuffer & buffer,
// 				  const Array<UInt> & dofs,
// 				  SynchronizationTag tag) {
//   AKANTU_DEBUG_IN();

//   switch(tag) {
//   case SynchronizationTag::_solver_solution: {
//     unpackDOFDataHelper(*solution, buffer, dofs);
//     break;
//   }
//   default: {
//   }
//   }

//   AKANTU_DEBUG_OUT();
// }
