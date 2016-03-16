/**
 * @file   sparse_matrix_aij.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Aug 18 16:31:23 2015
 *
 * @brief  Implementation of the AIJ sparse matrix
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

/* -------------------------------------------------------------------------- */
#include "sparse_matrix_aij.hh"
#include "dof_manager_default.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
SparseMatrixAIJ::SparseMatrixAIJ(DOFManagerDefault & dof_manager,
                                 const MatrixType & matrix_type, const ID & id)
    : SparseMatrix(dof_manager, matrix_type, id), dof_manager(dof_manager),
      irn(0, 1, id + ":irn"), jcn(0, 1, id + ":jcn"), a(0, 1, id + ":a"),
      profile_release(0), value_release(0) {}

/* -------------------------------------------------------------------------- */
SparseMatrixAIJ::SparseMatrixAIJ(const SparseMatrixAIJ & matrix, const ID & id)
    : SparseMatrix(matrix, id), dof_manager(matrix.dof_manager),
      irn(matrix.irn, true, id + ":irn"), jcn(matrix.jcn, true, id + ":jcn"),
      a(matrix.a, true, id + ":a"), profile_release(0), value_release(0) {}

/* -------------------------------------------------------------------------- */
SparseMatrixAIJ::~SparseMatrixAIJ() {}

// /* --------------------------------------------------------------------------
// */
// void SparseMatrixAIJ::buildProfile(const Mesh & mesh,
//                                    const DOFSynchronizer & dof_synchronizer,
//                                    UInt nb_degree_of_freedom) {
//   AKANTU_DEBUG_IN();

//   // if(irn_jcn_to_k) delete irn_jcn_to_k;
//   // irn_jcn_to_k = new std::map<std::pair<UInt, UInt>, UInt>;
//   clearProfile();

//   this->dof_synchronizer = &const_cast<DOFSynchronizer &>(dof_synchronizer);

//   coordinate_list_map::iterator irn_jcn_k_it;

//   Int * eq_nb_val = dof_synchronizer.getGlobalDOFEquationNumbers().storage();

//   Mesh::type_iterator it =
//       mesh.firstType(mesh.getSpatialDimension(), _not_ghost,
//       _ek_not_defined);
//   Mesh::type_iterator end =
//       mesh.lastType(mesh.getSpatialDimension(), _not_ghost, _ek_not_defined);
//   for (; it != end; ++it) {
//     UInt nb_element = mesh.getNbElement(*it);
//     UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
//     UInt size_mat = nb_nodes_per_element * nb_degree_of_freedom;

//     UInt * conn_val = mesh.getConnectivity(*it, _not_ghost).storage();
//     Int * local_eq_nb_val =
//         new Int[nb_degree_of_freedom * nb_nodes_per_element];

//     for (UInt e = 0; e < nb_element; ++e) {
//       Int * tmp_local_eq_nb_val = local_eq_nb_val;
//       for (UInt i = 0; i < nb_nodes_per_element; ++i) {
//         UInt n = conn_val[i];
//         for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
//           *tmp_local_eq_nb_val++ = eq_nb_val[n * nb_degree_of_freedom + d];
//         }
//         // memcpy(tmp_local_eq_nb_val, eq_nb_val + n * nb_degree_of_freedom,
//         // nb_degree_of_freedom * sizeof(Int));
//         // tmp_local_eq_nb_val += nb_degree_of_freedom;
//       }

//       for (UInt i = 0; i < size_mat; ++i) {
//         UInt c_irn = local_eq_nb_val[i];
//         if (c_irn < this->size) {
//           UInt j_start = (sparse_matrix_type == _symmetric) ? i : 0;
//           for (UInt j = j_start; j < size_mat; ++j) {
//             UInt c_jcn = local_eq_nb_val[j];
//             if (c_jcn < this->size) {
//               KeyCOO irn_jcn = key(c_irn, c_jcn);
//               irn_jcn_k_it = irn_jcn_k.find(irn_jcn);

//               if (irn_jcn_k_it == irn_jcn_k.end()) {
//                 irn_jcn_k[irn_jcn] = nb_non_zero;
//                 irn.push_back(c_irn + 1);
//                 jcn.push_back(c_jcn + 1);
//                 nb_non_zero++;
//               }
//             }
//           }
//         }
//       }
//       conn_val += nb_nodes_per_element;
//     }

//     delete[] local_eq_nb_val;
//   }

//   /// for pbc @todo correct it for parallel
//   if (StaticCommunicator::getStaticCommunicator().getNbProc() == 1) {
//     for (UInt i = 0; i < size; ++i) {
//       KeyCOO irn_jcn = key(i, i);
//       irn_jcn_k_it = irn_jcn_k.find(irn_jcn);
//       if (irn_jcn_k_it == irn_jcn_k.end()) {
//         irn_jcn_k[irn_jcn] = nb_non_zero;
//         irn.push_back(i + 1);
//         jcn.push_back(i + 1);
//         nb_non_zero++;
//       }
//     }
//   }

//   a.resize(nb_non_zero);

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::applyBoundary(Real block_val) {
  AKANTU_DEBUG_IN();

  const Array<bool> & blocked_dofs = this->dof_manager.getGlobalBlockedDOFs();
  Array<Int>::const_scalar_iterator irn_val = irn.begin();
  Array<Int>::const_scalar_iterator jcn_val = jcn.begin();
  Array<Real>::scalar_iterator a_val = a.begin();

  for (UInt i = 0; i < nb_non_zero; ++i) {
    UInt ni = this->dof_manager.globalToLocalEquationNumber(*irn_val - 1);
    UInt nj = this->dof_manager.globalToLocalEquationNumber(*jcn_val - 1);
    if (blocked_dofs(ni) || blocked_dofs(nj)) {
      if (*irn_val != *jcn_val) {
        *a_val = 0;
      } else {
        if (this->dof_manager.isLocalOrMasterDOF(ni)) {
          *a_val = block_val;
        } else {
          *a_val = 0;
        }
      }
    }
    ++irn_val;
    ++jcn_val;
    ++a_val;
  }

  this->value_release++;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::saveProfile(const std::string & filename) const {
  AKANTU_DEBUG_IN();

  std::ofstream outfile;
  outfile.open(filename.c_str());

  outfile << "%%MatrixMarket matrix coordinate pattern";

  if (this->matrix_type == _symmetric)
    outfile << " symmetric";
  else
    outfile << " general";
  outfile << std::endl;

  UInt m = this->size;
  outfile << m << " " << m << " " << this->nb_non_zero << std::endl;

  for (UInt i = 0; i < this->nb_non_zero; ++i) {
    outfile << this->irn.storage()[i] << " " << this->jcn.storage()[i] << " 1"
            << std::endl;
  }

  outfile.close();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::saveMatrix(const std::string & filename) const {
  AKANTU_DEBUG_IN();

  std::ofstream outfile;
  outfile.precision(std::numeric_limits<Real>::digits10);

  outfile.open(filename.c_str());

  outfile << "%%MatrixMarket matrix coordinate real";

  if (this->matrix_type == _symmetric)
    outfile << " symmetric";
  else
    outfile << " general";
  outfile << std::endl;

  outfile << this->size << " " << this->size << " " << this->nb_non_zero
          << std::endl;

  for (UInt i = 0; i < this->nb_non_zero; ++i) {
    outfile << this->irn(i) << " " << this->jcn(i) << " " << this->a(i)
            << std::endl;
  }

  outfile.close();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::matVecMul(const Array<Real> & x, Array<Real> & y,
                                Real alpha, Real beta) const {
  AKANTU_DEBUG_IN();

  y *= beta;

  Array<Int>::const_scalar_iterator i_it = this->irn.storage();
  Array<Int>::const_scalar_iterator j_it = this->jcn.storage();
  Array<Real>::const_scalar_iterator a_it = this->a.storage();
  Array<Real>::const_scalar_iterator x_it = x.storage();
  Array<Real>::scalar_iterator y_it = y.storage();

  for (UInt k = 0; k < this->nb_non_zero; ++k, ++i_it, ++j_it, ++a_it) {
    const Int & i = *i_it;
    const Int & j = *j_it;
    const Real & A = *a_it;

    y_it[i] += alpha * A * x_it[j];

    if ((this->matrix_type == _symmetric) && (i != j))
      y_it[j] += alpha * A * x_it[i];
  }

  //  if (dof_synchronizer)
  //  dof_synchronizer->reduceSynchronize<AddOperation>(vect);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::copyContent(const SparseMatrix & matrix) {
  AKANTU_DEBUG_IN();
  const SparseMatrixAIJ & mat = dynamic_cast<const SparseMatrixAIJ &>(matrix);
  AKANTU_DEBUG_ASSERT(nb_non_zero == mat.getNbNonZero(),
                      "The to matrix don't have the same profiles");
  memcpy(a.storage(), mat.getA().storage(), nb_non_zero * sizeof(Real));

  this->value_release++;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::add(const SparseMatrix & B, Real alpha) {
  Array<Real>::scalar_iterator a_it = this->a.begin();

  try {
    const SparseMatrixAIJ & B_aij = dynamic_cast<const SparseMatrixAIJ &>(B);
    Array<Real>::const_scalar_iterator b_it = B_aij.a.begin();
    for (UInt n = 0; n < this->nb_non_zero; ++n, ++a_it, ++b_it) {
       *a_it += alpha * *b_it;
    }
  } catch(...) {
    Array<Int>::const_scalar_iterator i_it = this->irn.begin();
    Array<Int>::const_scalar_iterator j_it = this->jcn.begin();
    for (UInt n = 0; n < this->nb_non_zero; ++n, ++a_it, ++i_it, ++j_it) {
      const Int & i = *i_it;
      const Int & j = *j_it;
      Real & A_ij = *a_it;
      A_ij += alpha * B(i - 1, j - 1);
    }
  }

  this->value_release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::clear() {
  a.set(0.);

  this->value_release++;
}

__END_AKANTU__
