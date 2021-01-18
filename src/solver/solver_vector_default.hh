/**
 * @file   solver_vector_default.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Jan 01 2019
 *
 * @brief A Documented file.
 *
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
#include "solver_vector.hh"
/* -------------------------------------------------------------------------- */
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLVER_VECTOR_DEFAULT_HH_
#define AKANTU_SOLVER_VECTOR_DEFAULT_HH_

namespace akantu {
class DOFManagerDefault;
} // namespace akantu

namespace akantu {

class SolverVectorArray : public SolverVector {
public:
  SolverVectorArray(DOFManagerDefault & dof_manager, const ID & id);
  SolverVectorArray(const SolverVectorArray & vector, const ID & id);

  ~SolverVectorArray() override = default;

  virtual Array<Real> & getVector() = 0;
  virtual const Array<Real> & getVector() const = 0;

  void printself(std::ostream & stream, int indent = 0) const override {
    std::string space(indent, AKANTU_INDENT);
    stream << space << "SolverVectorArray [" << std::endl;
    stream << space << " + id: " << id << std::endl;
    this->getVector().printself(stream, indent + 1);
    stream << space << "]" << std::endl;
  }
};

/* -------------------------------------------------------------------------- */
template <class Array_> class SolverVectorArrayTmpl : public SolverVectorArray {
public:
  SolverVectorArrayTmpl(DOFManagerDefault & dof_manager, Array_ & vector,
                        const ID & id = "solver_vector_default")
      : SolverVectorArray(dof_manager, id), dof_manager(dof_manager),
        vector(vector) {}

  template <class A = Array_,
            std::enable_if_t<not std::is_reference<A>::value> * = nullptr>
  SolverVectorArrayTmpl(DOFManagerDefault & dof_manager,
                        const ID & id = "solver_vector_default")
      : SolverVectorArray(dof_manager, id), dof_manager(dof_manager),
        vector(0, 1, id + ":vector") {}

  SolverVectorArrayTmpl(const SolverVectorArrayTmpl & vector,
                        const ID & id = "solver_vector_default")
      : SolverVectorArray(vector, id), dof_manager(vector.dof_manager),
        vector(vector.vector) {}

  operator const Array<Real> &() const override { return getVector(); };
  virtual operator Array<Real> &() { return getVector(); };

  SolverVector & operator+(const SolverVector & y) override;
  SolverVector & operator=(const SolverVector & y) override;

  void resize() override {
    static_assert(not std::is_const<std::remove_reference_t<Array_>>::value,
                  "Cannot resize a const Array");
    this->vector.resize(this->localSize(), 0.);
    ++this->release_;
  }

  void set(Real val) override {
      static_assert(not std::is_const<std::remove_reference_t<Array_>>::value,
                  "Cannot clear a const Array");
    this->vector.set(val);
    ++this->release_;
  }

public:
  Array<Real> & getVector() override { return vector; }
  const Array<Real> & getVector() const override { return vector; }

  Int size() const override;
  Int localSize() const override;

  virtual Array<Real> & getGlobalVector() { return this->vector; }
  virtual void setGlobalVector(const Array<Real> & solution) {
    this->vector.copy(solution);
  }

protected:
  DOFManagerDefault & dof_manager;
  Array_ vector;

  template <class A> friend class SolverVectorArrayTmpl;
};

/* -------------------------------------------------------------------------- */
using SolverVectorDefault = SolverVectorArrayTmpl<Array<Real>>;

/* -------------------------------------------------------------------------- */
template <class Array>
using SolverVectorDefaultWrap = SolverVectorArrayTmpl<Array &>;

template <class Array>
decltype(auto) make_solver_vector_default_wrap(DOFManagerDefault & dof_manager,
                                               Array & vector) {
  return SolverVectorDefaultWrap<Array>(dof_manager, vector);
}

} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "solver_vector_default_tmpl.hh"
/* -------------------------------------------------------------------------- */

#endif /* AKANTU_SOLVER_VECTOR_DEFAULT_HH_ */
