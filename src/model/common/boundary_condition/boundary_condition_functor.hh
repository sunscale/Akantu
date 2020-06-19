/**
 * @file   boundary_condition_functor.hh
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Definitions of the functors to apply boundary conditions
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
#include "aka_common.hh"
#include "fe_engine.hh"
#include "integration_point.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_BOUNDARY_CONDITION_FUNCTOR_HH__
#define __AKANTU_BOUNDARY_CONDITION_FUNCTOR_HH__

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
namespace BC {
  using Axis = ::akantu::SpatialDirection;

  /* ---------------------------------------------------------------------- */
  struct Functor {
    enum Type { _dirichlet, _neumann };
  };

  /* ---------------------------------------------------------------------- */
  namespace Dirichlet {

    class DirichletFunctor : public Functor {
    public:
      DirichletFunctor() = default;
      explicit DirichletFunctor(Axis ax) : axis(ax) {}

      virtual void operator()(__attribute__((unused)) UInt node,
                              __attribute__((unused)) Vector<bool> & flags,
                              __attribute__((unused)) Vector<Real> & primal,
                              __attribute__((unused))
                              const Vector<Real> & coord) const {
        AKANTU_TO_IMPLEMENT();
      }

    public:
      static const Type type = _dirichlet;

    protected:
      Axis axis{_x};
    };

    /* ---------------------------------------------------------------------- */
    class FlagOnly : public DirichletFunctor {
    public:
      explicit FlagOnly(Axis ax = _x) : DirichletFunctor(ax) {}

    public:
      inline void operator()(UInt node, Vector<bool> & flags,
                             Vector<Real> & primal,
                             const Vector<Real> & coord) const;
    };

    /* ---------------------------------------------------------------------- */
    // class FreeBoundary : public DirichletFunctor {
    // public:
    //   explicit FreeBoundary(Axis ax = _x) : DirichletFunctor(ax) {}

    // public:
    //   inline void operator()(UInt node, Vector<bool> & flags,
    //                          Vector<Real> & primal,
    //                          const Vector<Real> & coord) const;
    // };

    /* ---------------------------------------------------------------------- */
    class FixedValue : public DirichletFunctor {
    public:
      FixedValue(Real val, Axis ax = _x) : DirichletFunctor(ax), value(val) {}

    public:
      inline void operator()(UInt node, Vector<bool> & flags,
                             Vector<Real> & primal,
                             const Vector<Real> & coord) const;

    protected:
      Real value;
    };

    /* ---------------------------------------------------------------------- */
    class IncrementValue : public DirichletFunctor {
    public:
      IncrementValue(Real val, Axis ax = _x)
          : DirichletFunctor(ax), value(val) {}

    public:
      inline void operator()(UInt node, Vector<bool> & flags,
                             Vector<Real> & primal,
                             const Vector<Real> & coord) const;

      inline void setIncrement(Real val) { this->value = val; }

    protected:
      Real value;
    };

    /* ---------------------------------------------------------------------- */
    class Increment : public DirichletFunctor {
    public:
      explicit Increment(const Vector<Real> & val)
          : DirichletFunctor(_x), value(val) {}

    public:
      inline void operator()(UInt node, Vector<bool> & flags,
                             Vector<Real> & primal,
                             const Vector<Real> & coord) const;

      inline void setIncrement(const Vector<Real> & val) { this->value = val; }

    protected:
      Vector<Real> value;
    };
  } // namespace Dirichlet

  /* ------------------------------------------------------------------------ */
  /* Neumann                                                                  */
  /* ------------------------------------------------------------------------ */
  namespace Neumann {

    class NeumannFunctor : public Functor {

    protected:
      NeumannFunctor() = default;

    public:
      virtual void operator()(const IntegrationPoint & quad_point,
                              Vector<Real> & dual, const Vector<Real> & coord,
                              const Vector<Real> & normals) const = 0;

      virtual ~NeumannFunctor() = default;

    public:
      static const Type type = _neumann;
    };

    /* ---------------------------------------------------------------------- */
    class FromHigherDim : public NeumannFunctor {
    public:
      explicit FromHigherDim(const Matrix<Real> & mat) : bc_data(mat) {}
      ~FromHigherDim() override = default;

    public:
      inline void operator()(const IntegrationPoint & quad_point,
                             Vector<Real> & dual, const Vector<Real> & coord,
                             const Vector<Real> & normals) const override;

    protected:
      Matrix<Real> bc_data;
    };

    /* ---------------------------------------------------------------------- */
    class FromSameDim : public NeumannFunctor {
    public:
      explicit FromSameDim(const Vector<Real> & vec) : bc_data(vec) {}
      ~FromSameDim() override = default;

    public:
      inline void operator()(const IntegrationPoint & quad_point,
                             Vector<Real> & dual, const Vector<Real> & coord,
                             const Vector<Real> & normals) const override;

    protected:
      Vector<Real> bc_data;
    };

    /* ---------------------------------------------------------------------- */
    class FreeBoundary : public NeumannFunctor {
    public:
      inline void operator()(const IntegrationPoint & quad_point,
                             Vector<Real> & dual, const Vector<Real> & coord,
                             const Vector<Real> & normals) const override;
    };
  } // namespace Neumann
} // namespace BC
} // namespace akantu

#include "boundary_condition_functor_inline_impl.hh"

#endif /* __AKANTU_BOUNDARY_CONDITION_FUNCTOR_HH__ */
