/**
 * @file   integrator.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Dec 03 2017
 *
 * @brief  interface for integrator classes
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef __AKANTU_INTEGRATOR_HH__
#define __AKANTU_INTEGRATOR_HH__

/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {

class Integrator : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  Integrator(const Mesh & mesh, UInt spatial_dimension,
             const ID & id = "integrator", const MemoryID & memory_id = 0)
      : Memory(id, memory_id), mesh(mesh),
        _spatial_dimension(spatial_dimension),
        jacobians("jacobians", id, memory_id) {
    AKANTU_DEBUG_IN();

    AKANTU_DEBUG_OUT();
  };

  ~Integrator() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// empty method
  template <ElementType type>
  inline void precomputeJacobiansOnQuadraturePoints(__attribute__((unused))
                                                    GhostType ghost_type) {}

  /// empty method
  void integrateOnElement(const Array<Real> & /*f*/, Real * /*intf*/,
                          UInt /*nb_degree_of_freedom*/,
                          const Element & /*elem*/,
                          GhostType /*ghost_type*/) const {};

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
      ;
    stream << space << "Integrator [" << std::endl;
    jacobians.printself(stream, indent + 1);
    stream << space << "]" << std::endl;
  };

  /* ------------------------------------------------------------------------ */
public:
  virtual void onElementsAdded(const Array<Element> &) {}
  virtual void
  onElementsRemoved(const Array<Element> &,
                    const ElementTypeMapArray<UInt> & new_numbering) {
    jacobians.onElementsRemoved(new_numbering);
  }
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// access to the jacobians
  Array<Real> & getJacobians(const ElementType & type,
                             const GhostType & ghost_type = _not_ghost) {
    return jacobians(type, ghost_type);
  };

  /// access to the jacobians const
  const Array<Real> &
  getJacobians(const ElementType & type,
               const GhostType & ghost_type = _not_ghost) const {
    return jacobians(type, ghost_type);
  };

  AKANTU_GET_MACRO(Jacobians, jacobians, const ElementTypeMapArray<Real> &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// mesh associated to the integrator
  const Mesh & mesh;

  // spatial dimension of the elements to consider
  UInt _spatial_dimension;

  /// jacobians for all elements
  ElementTypeMapArray<Real> jacobians;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "integrator_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const Integrator & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* __AKANTU_INTEGRATOR_HH__ */
