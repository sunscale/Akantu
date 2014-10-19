/**
 * @file   boundary_condition_functor.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 *
 * @date   Wed Mar 18 11:30:00 2013
 *
 * @brief  XXX
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

#define DIRICHLET_SANITY_CHECK \
  AKANTU_DEBUG_ASSERT(coord.size() <= flags.size(),                     \
                      "The coordinates and flags vectors given to the boundary" \
                      << " condition functor have different sizes!");   \
  AKANTU_DEBUG_ASSERT(primal.size() <= coord.size(),                    \
                      "The primal vector and coordinates vector given"  \
                      << " to the boundary condition functor have different sizes!");

#define NEUMANN_SANITY_CHECK                                            \
  AKANTU_DEBUG_ASSERT(coord.size() <= normals.size(),                   \
                      "The coordinates and normals vectors given to the" \
                      << " boundary condition functor have different sizes!"); \
  AKANTU_DEBUG_ASSERT(dual.size() <= coord.size(),                      \
                      "The dual vector and coordinates vector given to" \
                      << " the boundary condition functor have different sizes!");

__BEGIN_AKANTU__

namespace BC {
  namespace Dirichlet {
    /* ---------------------------------------------------------------------- */
    inline void FlagOnly::operator()(UInt node,
                                     Vector<bool> & flags,
                                     Vector<Real> & primal,
                                     const Vector<Real> & coord) const {

      DIRICHLET_SANITY_CHECK;

      flags(axis) = true;
    }

    /* ---------------------------------------------------------------------- */
    inline void FreeBoundary::operator()(UInt node,
                                         Vector<bool> & flags,
                                         Vector<Real> & primal,
                                         const Vector<Real> & coord) const {

      DIRICHLET_SANITY_CHECK;

      flags(axis) = false;
    }

    /* ---------------------------------------------------------------------- */
    inline void FixedValue::operator()(UInt node,
                                       Vector<bool> & flags,
                                       Vector<Real> & primal,
                                       const Vector<Real> & coord) const {
      DIRICHLET_SANITY_CHECK;
      flags(axis) = true;
      primal(axis) = value;
    }

    /* ---------------------------------------------------------------------- */
    inline void IncrementValue::operator()(UInt node,
                                           Vector<bool> & flags,
                                           Vector<Real> & primal,
                                           const Vector<Real> & coord) const {
      DIRICHLET_SANITY_CHECK;
      flags(axis) = true;
      primal(axis) += value;
    }
  } // end namespace Dirichlet

  /* ------------------------------------------------------------------------ */
  /* Neumann                                                                  */
  /* ------------------------------------------------------------------------ */
  namespace Neumann {
    /* ---------------------------------------------------------------------- */
    inline void FreeBoundary::operator()(const QuadraturePoint & quad_point,
                                         Vector<Real> & dual,
                                         const Vector<Real> & coord,
                                         const Vector<Real> & normals) const {
      for(UInt i(0); i<dual.size(); ++i) {
        dual(i) = 0.0;
      }
    }

    /* -------------------------------------------------------------------------- */
    inline void FromHigherDim::operator()(const QuadraturePoint & quad_point,
                                          Vector<Real> & dual,
                                          const Vector<Real> & coord,
                                          const Vector<Real> & normals) const {
      dual.mul<false>(bc_data, normals);
    }

    /* -------------------------------------------------------------------------- */
    inline void FromSameDim::operator()(const QuadraturePoint & quad_point,
                                        Vector<Real> & dual,
                                        const Vector<Real> & coord,
                                        const Vector<Real> & normals) const {
      dual = bc_data;
    }
  } // end namespace Neumann
} // end namespace BC


__END_AKANTU__


