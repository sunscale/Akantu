/**
 * @file   sparse_matrix_aij.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Aug 17 21:22:57 2015
 *
 * @brief AIJ implementation of the SparseMatrix (this the format used by Mumps)
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
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SPARSE_MATRIX_AIJ_HH__
#define __AKANTU_SPARSE_MATRIX_AIJ_HH__

/* -------------------------------------------------------------------------- */
#ifndef __INTEL_COMPILER
namespace std {
  namespace tr1 {
    /// hashing function for pairs
    template <typename a, typename b> class hash<std::pair<a, b>> {
    public:
      hash() : ah(), bh() {}
      size_t operator()(const std::pair<a, b> & p) const {
        size_t seed = ah(p.first);
        return bh(p.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }

    private:
      const hash<a> ah;
      const hash<b> bh;
    };
  }
} // namespaces
#endif

__BEGIN_AKANTU__

class SparseMatrixAIJ : public SparseMatrix {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseMatrixAIJ(DOFManager & dof_manager, const MatrixType & matrix_type,
                  const ID & id = "sparse_matrix",
                  const MemoryID & memory_id = 0);

  SparseMatrixAIJ(const SparseMatrix & matrix, const ID & id = "sparse_matrix",
                  const MemoryID & memory_id = 0);

  virtual ~SparseMatrixAIJ();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// remove the existing profile
  inline void clearProfile();

  /// add a non-zero element
  virtual UInt addToProfile(UInt i, UInt j);

  /// set the matrix to 0
  virtual void clear();

  /// assemble a local matrix in the sparse one
  inline void addToMatrix(UInt i, UInt j, Real value);

  /// set the size of the matrix
  void resize(UInt size) { this->size = size; }

  /// modify the matrix to "remove" the blocked dof
  virtual void applyBoundary(const Array<bool> & boundary, Real block_val = 1.);

  /// save the profil in a file using the MatrixMarket file format
  virtual void saveProfile(const std::string & filename) const;

  /// save the matrix in a file using the MatrixMarket file format
  virtual void saveMatrix(const std::string & filename) const;

  /// copy assuming the profile are the same
  virtual void copyContent(const SparseMatrix & matrix);

  /// add matrix assuming the profile are the same
  virtual void add(const SparseMatrix & matrix, Real alpha);

  /* ------------------------------------------------------------------------ */
  /// accessor to A_{ij} - if (i, j) not present it returns 0
  inline Real operator()(UInt i, UInt j) const;

  /// accessor to A_{ij} - if (i, j) not present it fails, (i, j) should be
  /// first added to the profile
  inline Real & operator()(UInt i, UInt j);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(IRN, irn, const Array<Int> &);

  AKANTU_GET_MACRO(JCN, jcn, const Array<Int> &);

  AKANTU_GET_MACRO(A, a, const Array<Real> &);

protected:
  typedef std::pair<UInt, UInt> KeyCOO;
  typedef unordered_map<KeyCOO, UInt>::type coordinate_list_map;

  /// get the pair corresponding to (i, j)
  inline KeyCOO key(UInt i, UInt j) const {
    if (this->matrix_type == _symmetric && (i > j)) return std::make_pair(j, i);
    return std::make_pair(i, j);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// row indexes
  Array<Int> irn;

  /// column indexes
  Array<Int> jcn;

  /// values : A[k] = Matrix[irn[k]][jcn[k]]
  Array<Real> a;

  /*
   * map for (i, j) ->  k correspondence
   */
  coordinate_list_map irn_jcn_k;
};

__END_AKANTU__

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "sparse_matrix_aij_inline_impl.cc"

#endif /* __AKANTU_SPARSE_MATRIX_AIJ_HH__ */
