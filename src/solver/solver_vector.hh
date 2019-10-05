/**
 * @file   solver_vector.hh
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
#include "aka_array.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLVER_VECTOR_HH__
#define __AKANTU_SOLVER_VECTOR_HH__

namespace akantu {
class DOFManager;
}

namespace akantu {

class SolverVector {
public:
  SolverVector(DOFManager & dof_manager, const ID & id = "solver_vector")
      : id(id), _dof_manager(dof_manager) {}

  SolverVector(const SolverVector & vector, const ID & id = "solver_vector")
      : id(id), _dof_manager(vector._dof_manager) {}

  virtual ~SolverVector() = default;

  // resize the vector to the size of the problem
  virtual void resize() = 0;

  // clear the vector
  virtual void clear() = 0;

  virtual operator const Array<Real> &() const = 0;

  virtual Int size() const = 0;
  virtual Int localSize() const = 0;

  virtual SolverVector & operator+(const SolverVector & y) = 0;
  virtual SolverVector & operator=(const SolverVector & y) = 0;

  UInt & release() { return release_; }
  UInt release() const { return release_; }

  virtual void printself(std::ostream & stream, int indent = 0) const = 0;

protected:
  ID id;

  /// Underlying dof manager
  DOFManager & _dof_manager;

  UInt release_{0};
};

inline std::ostream & operator<<(std::ostream & stream, SolverVector & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* __AKANTU_SOLVER_VECTOR_HH__ */
