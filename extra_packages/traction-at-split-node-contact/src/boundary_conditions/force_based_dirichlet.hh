/**
 * @file   force_based_dirichlet.hh
 *
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  dirichlet boundary condition that tries
 * to keep the force at a given value
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
  ForceBasedDirichlet(SolidMechanicsModel & model, 
		      BC::Axis ax, 
		      Real target_f,
		      Real mass = 0.) :
    IncrementValue(0., ax), 
    model(model),
    mass(mass), 
    velocity(0.), 
    target_force(target_f),
    total_residual(0.)
  {}

  virtual ~ForceBasedDirichlet() {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void updateTotalResidual() {
    SubBoundarySet::iterator it  = this->subboundaries.begin();
    SubBoundarySet::iterator end = this->subboundaries.end();
    this->total_residual = 0.;
    for (; it != end; ++it) {
      this->total_residual += integrateResidual(*it, this->model, this->axis);
    }
  }

  virtual Real update() {
    AKANTU_DEBUG_IN();
    
    this->updateTotalResidual();
    Real total_force = this->target_force + this->total_residual;
    
    Real a = total_force / this->mass;
    Real dt = model.getTimeStep();  
    this->velocity += 0.5*dt*a;
    this->value     = this->velocity*dt + 0.5*dt*dt*a;  // increment position dx
    this->velocity += 0.5*dt*a;
    
    AKANTU_DEBUG_OUT();
    return this->total_residual;
  }

  Real applyYourself() {
    AKANTU_DEBUG_IN();
    Real reaction = this->update();

    SubBoundarySet::iterator it  = this->subboundaries.begin();
    SubBoundarySet::iterator end = this->subboundaries.end();
    for (; it != end; ++it) {
      this->model.applyBC(*this,*it);
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
