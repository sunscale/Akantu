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
#include "aka_iterators.hh"
#include "dof_manager_default.hh"
#include "dof_synchronizer.hh"
#include "terms_to_assemble.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SparseMatrixAIJ::SparseMatrixAIJ(DOFManagerDefault & dof_manager,
                                 const MatrixType & matrix_type, const ID & id)
    : SparseMatrix(dof_manager, matrix_type, id), dof_manager(dof_manager),
      irn(0, 1, id + ":irn"), jcn(0, 1, id + ":jcn"), a(0, 1, id + ":a") {}

/* -------------------------------------------------------------------------- */
SparseMatrixAIJ::SparseMatrixAIJ(const SparseMatrixAIJ & matrix, const ID & id)
    : SparseMatrix(matrix, id), dof_manager(matrix.dof_manager),
      irn(matrix.irn, true, id + ":irn"), jcn(matrix.jcn, true, id + ":jcn"),
      a(matrix.a, true, id + ":a") {}

/* -------------------------------------------------------------------------- */
SparseMatrixAIJ::~SparseMatrixAIJ() = default;

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::applyBoundary(Real block_val) {
  AKANTU_DEBUG_IN();

  // clang-format off
  const auto & blocked_dofs = this->dof_manager.getGlobalBlockedDOFs();

  for (auto && ij_a : zip(irn, jcn, a)) {
    UInt ni = this->dof_manager.globalToLocalEquationNumber(std::get<0>(ij_a) - 1);
    UInt nj = this->dof_manager.globalToLocalEquationNumber(std::get<1>(ij_a) - 1);
    if (blocked_dofs(ni) || blocked_dofs(nj)) {
      std::get<2>(ij_a) =
          std::get<0>(ij_a) != std::get<1>(ij_a)   ? 0.
        : this->dof_manager.isLocalOrMasterDOF(ni) ? block_val
        :                                            0.;
    }
  }

  this->value_release++;
  // clang-format on

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::saveProfile(const std::string & filename) const {
  AKANTU_DEBUG_IN();

  std::ofstream outfile;
  outfile.open(filename.c_str());

  UInt m = this->size_;
  outfile << "%%MatrixMarket matrix coordinate pattern";
  if (this->matrix_type == _symmetric)
    outfile << " symmetric";
  else
    outfile << " general";
  outfile << std::endl;
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

  // open and set the properties of the stream
  std::ofstream outfile;
  outfile.open(filename.c_str());
  outfile.precision(std::numeric_limits<Real>::digits10);

  // write header
  outfile << "%%MatrixMarket matrix coordinate real";
  if (this->matrix_type == _symmetric)
    outfile << " symmetric";
  else
    outfile << " general";
  outfile << std::endl;
  outfile << this->size_ << " " << this->size_ << " " << this->nb_non_zero
          << std::endl;

  // write content
  for (UInt i = 0; i < this->nb_non_zero; ++i) {
    outfile << this->irn(i) << " " << this->jcn(i) << " " << this->a(i)
            << std::endl;
  }

  // time to end
  outfile.close();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::matVecMul(const Array<Real> & x, Array<Real> & y,
                                Real alpha, Real beta) const {
  AKANTU_DEBUG_IN();

  y *= beta;

  auto i_it = this->irn.begin();
  auto j_it = this->jcn.begin();
  auto a_it = this->a.begin();
  auto a_end = this->a.end();
  auto x_it = x.begin_reinterpret(x.size() * x.getNbComponent());
  auto y_it = y.begin_reinterpret(x.size() * x.getNbComponent());

  for (; a_it != a_end; ++i_it, ++j_it, ++a_it) {
    Int i = this->dof_manager.globalToLocalEquationNumber(*i_it - 1);
    Int j = this->dof_manager.globalToLocalEquationNumber(*j_it - 1);
    const Real & A = *a_it;

    y_it[i] += alpha * A * x_it[j];

    if ((this->matrix_type == _symmetric) && (i != j))
      y_it[j] += alpha * A * x_it[i];
  }

  if (this->dof_manager.hasSynchronizer())
    this->dof_manager.getSynchronizer().reduceSynchronize<AddOperation>(y);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::copyContent(const SparseMatrix & matrix) {
  AKANTU_DEBUG_IN();
  const auto & mat = dynamic_cast<const SparseMatrixAIJ &>(matrix);
  AKANTU_DEBUG_ASSERT(nb_non_zero == mat.getNbNonZero(),
                      "The to matrix don't have the same profiles");
  memcpy(a.storage(), mat.getA().storage(), nb_non_zero * sizeof(Real));

  this->value_release++;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class MatrixType>
void SparseMatrixAIJ::addMeToTemplated(MatrixType & B, Real alpha) const {
  UInt i, j;
  Real A_ij;
  for (auto && tuple : zip(irn, jcn, a)) {
    std::tie(i, j, A_ij) = tuple;
    B.add(i - 1, j - 1, alpha * A_ij);
  }
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::addMeTo(SparseMatrix & B, Real alpha) const {
  if (auto * B_aij = dynamic_cast<SparseMatrixAIJ *>(&B)) {
    this->addMeToTemplated<SparseMatrixAIJ>(*B_aij, alpha);
  } else {
    //    this->addMeToTemplated<SparseMatrix>(*this, alpha);
  }
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::mul(Real alpha) {
  this->a *= alpha;
  this->value_release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixAIJ::clear() {
  a.set(0.);

  this->value_release++;
}

} // namespace akantu
