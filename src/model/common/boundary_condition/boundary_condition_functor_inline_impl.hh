/**
 * @file   boundary_condition_functor_inline_impl.hh
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Mon Feb 19 2018
 *
 * @brief  implementation of the BC::Functors
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#include "boundary_condition_functor.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_BOUNDARY_CONDITION_FUNCTOR_INLINE_IMPL_HH_
#define AKANTU_BOUNDARY_CONDITION_FUNCTOR_INLINE_IMPL_HH_

/* -------------------------------------------------------------------------- */
#define DIRICHLET_SANITY_CHECK                                                 \
  AKANTU_DEBUG_ASSERT(                                                         \
      coord.size() <= flags.size(),                                            \
      "The coordinates and flags vectors given to the boundary"                \
          << " condition functor have different sizes!");                      \
  AKANTU_DEBUG_ASSERT(                                                         \
      primal.size() <= coord.size(),                                           \
      "The primal vector and coordinates vector given"                         \
          << " to the boundary condition functor have different sizes!");

#define NEUMANN_SANITY_CHECK                                                   \
  AKANTU_DEBUG_ASSERT(                                                         \
      coord.size() <= normals.size(),                                          \
      "The coordinates and normals vectors given to the"                       \
          << " boundary condition functor have different sizes!");             \
  AKANTU_DEBUG_ASSERT(                                                         \
      dual.size() <= coord.size(),                                             \
      "The dual vector and coordinates vector given to"                        \
          << " the boundary condition functor have different sizes!");

namespace akantu {
namespace BC {
  /* ---------------------------------------------------------------------- */
  namespace Dirichlet {
    inline void FlagOnly::
    operator()(__attribute__((unused)) UInt node, Vector<bool> & flags,
               __attribute__((unused)) Vector<Real> & primal,
               __attribute__((unused)) const Vector<Real> & coord) const {

      DIRICHLET_SANITY_CHECK;

      flags(this->axis) = true;
    }

    /* ---------------------------------------------------------------------- */
    // inline void FreeBoundary::
    // operator()(__attribute__((unused)) UInt node, Vector<bool> & flags,
    //            __attribute__((unused)) Vector<Real> & primal,
    //            __attribute__((unused)) const Vector<Real> & coord) const {

    //   DIRICHLET_SANITY_CHECK;

    //   flags(this->axis) = false;
    // }

    /* ---------------------------------------------------------------------- */
    inline void FixedValue::operator()(__attribute__((unused)) UInt node,
                                       Vector<bool> & flags,
                                       Vector<Real> & primal,
                                       __attribute__((unused))
                                       const Vector<Real> & coord) const {
      DIRICHLET_SANITY_CHECK;
      flags(this->axis) = true;
      primal(this->axis) = value;
    }

    /* ---------------------------------------------------------------------- */
    inline void IncrementValue::operator()(__attribute__((unused)) UInt node,
                                           Vector<bool> & flags,
                                           Vector<Real> & primal,
                                           __attribute__((unused))
                                           const Vector<Real> & coord) const {
      DIRICHLET_SANITY_CHECK;
      flags(this->axis) = true;
      primal(this->axis) += value;
    }

    /* ---------------------------------------------------------------------- */
    inline void Increment::operator()(__attribute__((unused)) UInt node,
                                      Vector<bool> & flags,
                                      Vector<Real> & primal,
                                      __attribute__((unused))
                                      const Vector<Real> & coord) const {
      DIRICHLET_SANITY_CHECK;
      flags.set(true);
      primal += value;
    }
  } // namespace Dirichlet
  /* ------------------------------------------------------------------------ */
  /* Neumann */
  /* ------------------------------------------------------------------------ */

  namespace Neumann {
    inline void FreeBoundary::
    operator()(__attribute__((unused)) const IntegrationPoint & quad_point,
               Vector<Real> & dual,
               __attribute__((unused)) const Vector<Real> & coord,
               __attribute__((unused)) const Vector<Real> & normals) const {
      for (UInt i(0); i < dual.size(); ++i) {
        dual(i) = 0.0;
      }
    }

    /* ---------------------------------------------------------------------- */
    inline void FromHigherDim::operator()(__attribute__((unused))
                                          const IntegrationPoint & quad_point,
                                          Vector<Real> & dual,
                                          __attribute__((unused))
                                          const Vector<Real> & coord,
                                          const Vector<Real> & normals) const {
      dual.mul<false>(this->bc_data, normals);
    }

    /* ---------------------------------------------------------------------- */
    inline void FromSameDim::
    operator()(__attribute__((unused)) const IntegrationPoint & quad_point,
               Vector<Real> & dual,
               __attribute__((unused)) const Vector<Real> & coord,
               __attribute__((unused)) const Vector<Real> & normals) const {
      dual = this->bc_data;
    }
  } // namespace Neumann
} // namespace BC
} // namespace akantu

#endif /* AKANTU_BOUNDARY_CONDITION_FUNCTOR_INLINE_IMPL_HH_ */
