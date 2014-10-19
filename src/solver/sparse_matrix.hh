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
#include "aka_common.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */
#ifndef __INTEL_COMPILER
namespace std {
  namespace tr1 {
    template<typename a, typename b>
    struct hash< std::pair<a, b> > {
    private:
      const hash<a> ah;
      const hash<b> bh;
    public:
      hash() : ah(), bh() {}
      size_t operator()(const std::pair<a, b> &p) const {
	size_t seed = ah(p.first);
	return bh(p.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);
      }
    };
  }
} // namespaces
#endif

__BEGIN_AKANTU__

class DOFSynchronizer;

class SparseMatrix : private Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseMatrix(UInt size,
	       const SparseMatrixType & sparse_matrix_type,
	       const ID & id = "sparse_matrix",
	       const MemoryID & memory_id = 0);

  SparseMatrix(const SparseMatrix & matrix,
	       const ID & id = "sparse_matrix",
	       const MemoryID & memory_id = 0);

  virtual ~SparseMatrix();

  typedef std::pair<UInt, UInt> KeyCOO;
  typedef unordered_map<KeyCOO, UInt>::type coordinate_list_map;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// remove the existing profile
  inline void clearProfile();
  
  /// add a non-zero element
  virtual UInt addToProfile(UInt i, UInt j);

  /// set the matrix to 0
  inline void clear();

  /// assemble a local matrix in the sparse one
  inline void addToMatrix(UInt i, UInt j, Real value);

  /// set the size of the matrix
  void resize(UInt size)
  { this->size = size; }
  
  void buildProfile(const Mesh & mesh, const DOFSynchronizer & dof_synchronizer, UInt nb_degree_of_freedom);

  /// modify the matrix to "remove" the blocked dof
  virtual void applyBoundary(const Array<bool> & boundary, Real block_val = 1.);

//  /// modify the matrix to "remove" the blocked dof
//  void applyBoundaryNormal(Array<bool> & boundary_normal, Array<Real> & EulerAngles, Array<Real> & rhs, const Array<Real> & matrix, Array<Real> & rhs_rotated);

  /// modify the matrix to "remove" the blocked dof
  virtual void removeBoundary(const Array<bool> & boundary);

  /// restore the profile that was before removing the boundaries
  virtual void restoreProfile();

  /// save the profil in a file using the MatrixMarket file format
  virtual void saveProfile(const std::string & filename) const;

  /// save the matrix in a file using the MatrixMarket file format
  virtual void saveMatrix(const std::string & filename) const;

  /// copy assuming the profile are the same
  virtual void copyContent(const SparseMatrix & matrix);
  
  /// copy profile
//  void copyProfile(const SparseMatrix & matrix);

  /// add matrix assuming the profile are the same
  virtual void add(const SparseMatrix & matrix, Real alpha);

  /// diagonal lumping
  virtual void lump(Array<Real> & lumped);

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  inline KeyCOO key(UInt i, UInt j) const {
    if(sparse_matrix_type == _symmetric && (i > j))
      return std::make_pair(j, i);

    return std::make_pair(i, j);
  }


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the values at potition i, j
  inline Real operator()(UInt i, UInt j) const;
  inline Real & operator()(UInt i, UInt j);

  AKANTU_GET_MACRO(IRN, irn, const Array<Int> &);

  AKANTU_GET_MACRO(JCN, jcn, const Array<Int> &);

  AKANTU_GET_MACRO(A, a, const Array<Real> &);

  AKANTU_GET_MACRO(NbNonZero, nb_non_zero, UInt);

  AKANTU_GET_MACRO(Size, size, UInt);

  AKANTU_GET_MACRO(SparseMatrixType, sparse_matrix_type, const SparseMatrixType &);

  const DOFSynchronizer & getDOFSynchronizer() const {
    AKANTU_DEBUG_ASSERT(dof_synchronizer != NULL,
			"DOFSynchronizer not initialized in the SparseMatrix!");
    return *dof_synchronizer;
  }

private:
  AKANTU_GET_MACRO(DOFSynchronizerPointer, dof_synchronizer, DOFSynchronizer *);

  friend Array<Real> & operator*=(Array<Real> & vect, const SparseMatrix & mat);
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id  of the SparseMatrix
  ID id;

  /// sparce matrix type
  SparseMatrixType sparse_matrix_type;

  /// Mesh corresponding to the profile
  //  const Mesh * mesh;

  /// Size of the matrix
  UInt size;

  /// number of processors
  UInt nb_proc;

  /// number of non zero element
  UInt nb_non_zero;

  /// row indexes
  Array<Int> irn;

  /// column indexes
  Array<Int> jcn;

  /// values : A[k] = Matrix[irn[k]][jcn[k]]
  Array<Real> a;


  /// saved row indexes
  Array<Int> * irn_save;

  /// saved column indexes
  Array<Int> * jcn_save;

  /// saved size
  UInt size_save;

  /// information to know where to assemble an element in a global sparse matrix
  //  ElementTypeMapArray<UInt> element_to_sparse_profile;

  /* map for  (i,j) ->  k correspondence \warning std::map are slow
   *  \todo improve  with hash_map (non standard in stl) or unordered_map (boost or C++0x)
   */
  coordinate_list_map irn_jcn_k;

  DOFSynchronizer * dof_synchronizer;
  //  std::map<std::pair<UInt, UInt>, UInt> * irn_jcn_to_k;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "sparse_matrix_inline_impl.cc"
#endif

// /// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const SparseMatrix & _this)
// {
//   _this.printself(stream);
//   return stream;
// }

Array<Real> & operator*=(Array<Real> & vect, const SparseMatrix & mat);


__END_AKANTU__

#endif /* __AKANTU_SPARSE_MATRIX_HH__ */
