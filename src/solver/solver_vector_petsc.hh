/**
 * @file   solver_vector_petsc.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Jan 01 2019
 *
 * @brief A Documented file.
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
#include "dof_manager_petsc.hh"
#include "solver_vector.hh"
/* -------------------------------------------------------------------------- */
#include <petscvec.h>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLVER_VECTOR_PETSC_HH__
#define __AKANTU_SOLVER_VECTOR_PETSC_HH__

namespace akantu {
class DOFManagerPETSc;
} // namespace akantu

namespace akantu {

/* -------------------------------------------------------------------------- */
namespace internal {
  /* ------------------------------------------------------------------------ */
  class PETScVector {
  public:
    virtual ~PETScVector() = default;

    operator Vec &() { return x; }
    operator const Vec &() const { return x; }

    Int size() const {
      PetscInt n;
      PETSc_call(VecGetSize, x, &n);
      return n;
    }
    Int local_size() const {
      PetscInt n;
      PETSc_call(VecGetLocalSize, x, &n);
      return n;
    }

    AKANTU_GET_MACRO_NOT_CONST(Vec, x, auto &);
    AKANTU_GET_MACRO(Vec, x, const auto &);

  protected:
    Vec x{nullptr};
  };

} // namespace internal

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
class SolverVectorPETSc : public SolverVector, public internal::PETScVector {
public:
  SolverVectorPETSc(DOFManagerPETSc & dof_manager,
                    const ID & id = "solver_vector_petsc");

  SolverVectorPETSc(const SolverVectorPETSc & vector,
                    const ID & id = "solver_vector_petsc");

  SolverVectorPETSc(Vec vec, DOFManagerPETSc & dof_manager,
                    const ID & id = "solver_vector_petsc");

  ~SolverVectorPETSc() override;

  // resize the vector to the size of the problem
  void resize() override;
  void clear() override;

  operator const Array<Real> &() const override;

  SolverVector & operator+(const SolverVector & y) override;
  SolverVector & operator=(const SolverVector & y) override;
  SolverVectorPETSc & operator=(const SolverVectorPETSc & y);

  /// get values using processors global indexes
  void getValues(const Array<Int> & idx, Array<Real> & values) const;

  /// get values using processors local indexes
  void getValuesLocal(const Array<Int> & idx, Array<Real> & values) const;

  /// adding values to the vector using the global indices
  void addValues(const Array<Int> & gidx, const Array<Real> & values,
                 Real scale_factor = 1.);

  /// adding values to the vector using the local indices
  void addValuesLocal(const Array<Int> & lidx, const Array<Real> & values,
                      Real scale_factor = 1.);

  Int size() const override { return internal::PETScVector::size(); }
  Int localSize() const override { return internal::PETScVector::local_size(); }

  void printself(std::ostream & stream, int indent = 0) const override;

protected:
  void applyModifications();
  void updateGhost();

protected:
  DOFManagerPETSc & dof_manager;

  // used for the conversion operator
  Array<Real> cache;
};

/* -------------------------------------------------------------------------- */
namespace internal {
  /* ------------------------------------------------------------------------ */
  template <class Array> class PETScWrapedVector : public PETScVector {
  public:
    PETScWrapedVector(Array && array) : array(array) {
      PETSc_call(VecCreateSeqWithArray, PETSC_COMM_SELF, 1, array.size(),
                 array.storage(), &x);
    }

    ~PETScWrapedVector() { PETSc_call(VecDestroy, &x); }

  private:
    Array array;
  };

  /* ------------------------------------------------------------------------ */
  template <bool read_only> class PETScLocalVector : public PETScVector {
  public:
    PETScLocalVector(const Vec & g) : g(g) {
      PETSc_call(VecGetLocalVectorRead, g, x);
    }
    PETScLocalVector(const SolverVectorPETSc & g)
        : PETScLocalVector(g.getVec()) {}
    ~PETScLocalVector() {
      PETSc_call(VecRestoreLocalVectorRead, g, x);
      PETSc_call(VecDestroy, &x);
    }

  private:
    const Vec & g;
  };

  template <> class PETScLocalVector<false> : public PETScVector {
  public:
    PETScLocalVector(Vec & g) : g(g) {
      PETSc_call(VecGetLocalVectorRead, g, x);
    }
    PETScLocalVector(SolverVectorPETSc & g) : PETScLocalVector(g.getVec()) {}
    ~PETScLocalVector() {
      PETSc_call(VecRestoreLocalVectorRead, g, x);
      PETSc_call(VecDestroy, &x);
    }

  private:
    Vec & g;
  };

  /* ------------------------------------------------------------------------ */
  template <class Array>
  decltype(auto) make_petsc_wraped_vector(Array && array) {
    return PETScWrapedVector<Array>(std::forward<Array>(array));
  }

  template <
      typename V,
      std::enable_if_t<std::is_same<Vec, std::decay_t<V>>::value> * = nullptr>
  decltype(auto) make_petsc_local_vector(V && vec) {
    constexpr auto read_only = std::is_const<std::remove_reference_t<V>>::value;
    return PETScLocalVector<read_only>(vec);
  }

  template <typename V, std::enable_if_t<std::is_base_of<
                            SolverVector, std::decay_t<V>>::value> * = nullptr>
  decltype(auto) make_petsc_local_vector(V && vec) {
    constexpr auto read_only = std::is_const<std::remove_reference_t<V>>::value;
    return PETScLocalVector<read_only>(
        dynamic_cast<std::conditional_t<read_only, const SolverVectorPETSc,
                                        SolverVectorPETSc> &>(vec));
  }

} // namespace internal

} // namespace akantu

#endif /* __AKANTU_SOLVER_VECTOR_PETSC_HH__ */
