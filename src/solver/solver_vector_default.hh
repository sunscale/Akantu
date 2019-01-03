/**
 * @file   solver_vector_default.hh
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
#include "solver_vector.hh"
/* -------------------------------------------------------------------------- */
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLVER_VECTOR_DEFAULT_HH__
#define __AKANTU_SOLVER_VECTOR_DEFAULT_HH__

namespace akantu {
class DOFManagerDefault;
} // namespace akantu

namespace akantu {

class SolverVectorArray : public SolverVector {
public:
  SolverVectorArray(DOFManagerDefault & dof_manager,
                    const ID & id = "solver_vector_default");

  SolverVectorArray(const SolverVectorArray & vector,
                    const ID & id = "solver_vector_default");

public:
  virtual Array<Real> & getVector() = 0;
  virtual const Array<Real> & getVector() const = 0;

protected:
  DOFManagerDefault & dof_manager;
};


/* -------------------------------------------------------------------------- */
class SolverVectorDefault : public SolverVectorArray {
public:
  SolverVectorDefault(DOFManagerDefault & dof_manager,
                      const ID & id = "solver_vector_default");

  SolverVectorDefault(const SolverVectorDefault & vector,
                      const ID & id = "solver_vector_default");

  // resize the vector to the size of the problem
  void resize() override;
  void clear() override;

public:
  Array<Real> & getVector() override { return vector; }
  const Array<Real> & getVector() const override { return vector; }

  virtual Array<Real> & getGlobalVector();
  virtual void setGlobalVector(const Array<Real> & global_vector);

protected:
  Array<Real> vector;
};

/* -------------------------------------------------------------------------- */
template <class Array> class SolverVectorDefaultWrap : public SolverVectorArray {
public:
  SolverVectorDefaultWrap(DOFManagerDefault & dof_manager, Array & vector);

  SolverVectorDefaultWrap(const SolverVectorDefaultWrap & vector,
                          const ID & id = "solver_vector_default");

  SolverVectorDefaultWrap(SolverVectorDefaultWrap && vector) = default;

  // resize the vector to the size of the problem
  void resize() override;

  // clear the vector
  void clear() override;

public:
  Array & getVector() override { return vector; }
  const Array & getVector() const override { return vector; }

protected:
  Array & vector;
};

template <class Array>
decltype(auto) make_solver_vector_default_wrap(DOFManagerDefault & dof_manager,
                                               Array & vector) {
  return SolverVectorDefaultWrap<Array>(dof_manager,
                                        vector);
}

} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "solver_vector_default_tmpl.hh"
/* -------------------------------------------------------------------------- */

#endif /* __AKANTU_SOLVER_VECTOR_DEFAULT_HH__ */
