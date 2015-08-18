/**
 * @file   sparse_matrix.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  sparse matrix storage class (distributed assembled matrix)
 * This is a COO format (Coordinate List)
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

#ifndef __AKANTU_SPARSE_MATRIX_HH__
#define __AKANTU_SPARSE_MATRIX_HH__

/* -------------------------------------------------------------------------- */
#include "dof_manager.hh"
#include "mesh.hh"

__BEGIN_AKANTU__

class SparseMatrix : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseMatrix(DOFManager & dof_manager,
               const MatrixType & matrix_type,
	       const ID & id = "sparse_matrix",
	       const MemoryID & memory_id = 0);

  SparseMatrix(const SparseMatrix & matrix,
	       const ID & id = "sparse_matrix",
	       const MemoryID & memory_id = 0);

  virtual ~SparseMatrix();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// remove the existing profile
  virtual void clearProfile();

  /// set the matrix to 0
  virtual void clear();

  /// add a non-zero element to the profile
  virtual inline UInt addToProfile(UInt i, UInt j) = 0;

  /// assemble a local matrix in the sparse one
  virtual inline void addToMatrix(UInt i, UInt j, Real value) = 0;

  /// save the profil in a file using the MatrixMarket file format
  virtual void saveProfile(const std::string & filename) const;

  /// save the matrix in a file using the MatrixMarket file format
  virtual void saveMatrix(const std::string & filename) const;

  /// copy assuming the profile are the same
  virtual void copyContent(const SparseMatrix & matrix);

  /// add matrix assuming the profile are the same
  virtual void add(const SparseMatrix & matrix, Real alpha);

  /// modify the matrix to "remove" the blocked dof
  virtual void applyBoundary(const Array<bool> & boundary, Real block_val = 1.);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the values at potition i, j
  virtual inline const Real & operator()(UInt i, UInt j) const { AKANTU_DEBUG_TO_IMPLEMENT(); }
  /// return the values at potition i, j
  virtual inline Real & operator()(UInt i, UInt j) { AKANTU_DEBUG_TO_IMPLEMENT(); }

  AKANTU_GET_MACRO(NbNonZero, nb_non_zero, UInt);

  AKANTU_GET_MACRO(Size, size, UInt);

  AKANTU_GET_MACRO(MatrixType, matrix_type, const MatrixType &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Underlying dof manager
  DOFManager & dof_manager;

  /// sparce matrix type
  MatrixType matrix_type;

  /// Size of the matrix
  UInt size;

  /// number of processors
  UInt nb_proc;

  /// number of non zero element
  UInt nb_non_zero;
};

Array<Real> & operator*=(Array<Real> & vect, const SparseMatrix & mat);

__END_AKANTU__

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "sparse_matrix_inline_impl.cc"


#endif /* __AKANTU_SPARSE_MATRIX_HH__ */
