/**
 * @file   petsc_matrix.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Aug 21 2015
 *
 * @brief  Interface for PETSc matrices
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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

#ifndef __AKANTU_PETSC_MATRIX_HH__
#define __AKANTU_PETSC_MATRIX_HH__

/* -------------------------------------------------------------------------- */
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include <petscao.h>
#include <petscmat.h>
/* -------------------------------------------------------------------------- */

namespace akantu {
class DOFManagerPETSc;
}

namespace akantu {

class SparseMatrixPETSc : public SparseMatrix {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseMatrixPETSc(DOFManagerPETSc & dof_manager,
                    const MatrixType & matrix_type,
                    const ID & id = "sparse_matrix",
                    const MemoryID & memory_id = 0);

  SparseMatrixPETSc(const SparseMatrix & matrix,
                    const ID & id = "sparse_matrix_petsc",
                    const MemoryID & memory_id = 0);

  virtual ~SparseMatrixPETSc();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// set the matrix to 0
  virtual void clear();

  /// modify the matrix to "remove" the blocked dof
  virtual void applyBoundary(const Array<bool> & boundary, Real block_val = 1.);

  /// save the matrix in a ASCII file format
  virtual void saveMatrix(const std::string & filename) const;

  /// add a sparse matrix assuming the profile are the same
  virtual void add(const SparseMatrix & matrix, Real alpha);
  /// add a petsc matrix assuming the profile are the same
  virtual void add(const SparseMatrixPETSc & matrix, Real alpha);

  Real operator()(UInt i, UInt j) const;

protected:
  typedef std::pair<UInt, UInt> KeyCOO;
  inline KeyCOO key(UInt i, UInt j) const { return std::make_pair(i, j); }

private:
  virtual void destroyInternalData();

  /// set the size of the PETSc matrix
  void setSize();
  void createGlobalAkantuToPETScMap(Int * local_master_eq_nbs_ptr);
  void createLocalAkantuToPETScMap();

  /// start to assemble the matrix
  void beginAssembly();
  /// finishes to assemble the matrix
  void endAssembly();

  /// perform the assembly stuff from petsc
  void performAssembly();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(PETScMat, mat, const Mat &);
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // DOFManagerPETSc that contains the numbering for petsc
  DOFManagerPETSc & dof_manager;

  /// store the PETSc matrix
  Mat mat;

  AO ao;

  /// size of the diagonal part of the matrix partition
  Int local_size;

  /// number of nonzeros in every row of the diagonal part
  Array<Int> d_nnz;

  /// number of nonzeros in every row of the off-diagonal part
  Array<Int> o_nnz;

  /// the global index of the first local row
  Int first_global_index;

  /// bool to indicate if the matrix data has been initialized by calling
  /// MatCreate
  bool is_petsc_matrix_initialized;
};

} // akantu

#endif /* __AKANTU_PETSC_MATRIX_HH__ */
