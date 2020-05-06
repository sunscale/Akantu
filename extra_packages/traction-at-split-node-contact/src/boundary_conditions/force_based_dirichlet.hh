/**
 * @file   force_based_dirichlet.hh
 *
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  dirichlet boundary condition that tries
 * to keep the force at a given value
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef __AST_FORCE_BASED_DIRICHLET_HH__
#define __AST_FORCE_BASED_DIRICHLET_HH__

// akantu
#include "aka_common.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
class ForceBasedDirichlet : public BC::Dirichlet::IncrementValue {

protected:
  typedef const Array<Real> * RealArrayPtr;
  typedef const Array<Int> * IntArrayPtr;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ForceBasedDirichlet(SolidMechanicsModel & model, BC::Axis ax, Real target_f,
                      Real mass = 0.)
      : IncrementValue(0., ax), model(model), mass(mass), velocity(0.),
        target_force(target_f), total_residual(0.) {}

  virtual ~ForceBasedDirichlet() {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void updateTotalResidual() {
    this->total_residual = 0.;
    for (auto && subboundary : this->subboundaries) {
      this->total_residual +=
          integrateResidual(subboundary, this->model, this->axis);
    }
  }

  virtual Real update() {
    AKANTU_DEBUG_IN();

    this->updateTotalResidual();
    Real total_force = this->target_force + this->total_residual;

    Real a = total_force / this->mass;
    Real dt = model.getTimeStep();
    this->velocity += 0.5 * dt * a;
    this->value =
        this->velocity * dt + 0.5 * dt * dt * a; // increment position dx
    this->velocity += 0.5 * dt * a;

    AKANTU_DEBUG_OUT();
    return this->total_residual;
  }

  Real applyYourself() {
    AKANTU_DEBUG_IN();
    Real reaction = this->update();

    for (auto && subboundary : this->subboundaries) {
      this->model.applyBC(*this, subboundary);
    }

    AKANTU_DEBUG_OUT();
    return reaction;
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_SET_MACRO(Mass, mass, Real);
  AKANTU_SET_MACRO(TargetForce, target_force, Real);

  void insertSubBoundary(const std::string & sb_name) {
    this->subboundaries.insert(sb_name);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  typedef std::set<std::string> SubBoundarySet;

protected:
  SolidMechanicsModel & model;
  SubBoundarySet subboundaries;

  Real mass;
  Real velocity;
  Real target_force;
  Real total_residual;
};

} // namespace akantu

#endif /* __AST_FORCE_BASED_DIRICHLET_HH__ */
