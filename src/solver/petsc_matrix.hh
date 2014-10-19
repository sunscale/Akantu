/**
 * @file   petsc_matrix.hh
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Jul 21 14:49:49 2014
 *
 * @brief  Interface for PETSc matrices
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

#ifndef __AKANTU_PETSC_MATRIX_HH__
#define __AKANTU_PETSC_MATRIX_HH__

/* -------------------------------------------------------------------------- */
#include "sparse_matrix.hh"

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

class PETScWrapper;

class PETScMatrix : public SparseMatrix {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PETScMatrix(UInt size,
	      const SparseMatrixType & sparse_matrix_type,
	      const ID & id = "petsc_matrix",
	      const MemoryID & memory_id = 0);

  PETScMatrix(const SparseMatrix & matrix,
	      const ID & id = "petsc_matrix",
	      const MemoryID & memory_id = 0);

  virtual ~PETScMatrix();


private:
  /// init internal PETSc matrix
  void init();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// resize PETSc matrix
  void resize(const DOFSynchronizer & dof_synchronizer);

  // /// remove the existing profile
  // inline void clearProfile();

  // /// add a non-zero element
  // virtual UInt addToProfile(UInt i, UInt j);

  // /// set the matrix to 0
  // inline void clear();

  // /// assemble a local matrix in the sparse one
  // inline void addToMatrix(UInt i, UInt j, Real value);

  /// fill the profil of the matrix
  virtual void buildProfile(const Mesh & mesh, const DOFSynchronizer & dof_synchronizer, UInt nb_degree_of_freedom);

  // /// modify the matrix to "remove" the blocked dof
  // virtual void applyBoundary(const Array<bool> & boundary, Real block_val = 1.);

  /// modify the matrix to "remove" the blocked dof
  virtual void applyBoundaryNormal(Array<bool> & boundary_normal, Array<Real> & EulerAngles, Array<Real> & rhs, const Array<Real> & matrix, Array<Real> & rhs_rotated) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// modify the matrix to "remove" the blocked dof
  virtual void removeBoundary(const Array<bool> & boundary) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
  
  /// perform assembly so that matrix is ready for use
  void performAssembly();

  // /// restore the profile that was before removing the boundaries
  // virtual void restoreProfile();

  // /// save the profil in a file using the MatrixMarket file format
  // virtual void saveProfile(const std::string & filename) const;

  /// save the matrix in a file using the MatrixMarket file format
  virtual void saveMatrix(const std::string & filename) const;

  // /// copy assuming the profile are the same
  // virtual void copyContent(const SparseMatrix & matrix);

  // /// copy profile
  // //  void copyProfile(const SparseMatrix & matrix);

  /// add matrix assuming the profile are the same
  virtual void add(const SparseMatrix & matrix, Real alpha);

  // /// diagonal lumping
  // virtual void lump(Array<Real> & lumped);

  // /// function to print the contain of the class
  // //virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  /// create a mapping from the global akantu numbering to the PETSc numbering
  void createGlobalAkantuToPETScMap(Int* local_master_eq_nbs_ptr);

  /// create a mapping from the local akantu numbering to the PETSc numbering
  void createLocalAkantuToPETScMap(const DOFSynchronizer & dof_synchronizer);


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */


private:
  /// store the PETSc structures
  PETScWrapper * petsc_wrapper;

  /// size of the diagonal part of the matrix partition
  Int local_size;

  /// number of nonzeros in every row of the diagonal part
  Array<Int> d_nnz;

  /// number of nonzeros in every row of the off-diagonal part
  Array<Int> o_nnz;
};

__END_AKANTU__

#endif /* __AKANTU_PETSCI_MATRIX_HH__ */
