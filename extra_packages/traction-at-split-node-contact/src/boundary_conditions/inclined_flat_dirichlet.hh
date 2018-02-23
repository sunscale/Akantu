/**
 * @file   inclined_flat_dirichlet.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  inclined dirichlet
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AST_INCLINED_FLAT_DIRICHLET_HH__
#define __AST_INCLINED_FLAT_DIRICHLET_HH__

// akantu
#include "aka_common.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
class InclinedFlatDirichlet : public BC::Dirichlet::DirichletFunctor {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  InclinedFlatDirichlet(Real val, BC::Axis ax, BC::Axis incl_ax,
                        Real center_coord, Real tang)
      : DirichletFunctor(ax), value(val), incl_ax(incl_ax),
        center_coord(center_coord), tang(tang){};

  virtual ~InclinedFlatDirichlet() {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline void operator()(UInt node, Vector<bool> & flags, Vector<Real> & primal,
                         const Vector<Real> & coord) const {
    AKANTU_DEBUG_IN();

    Real dist = coord(incl_ax) - this->center_coord;
    flags(axis) = true;
    primal(axis) = this->value + this->tang * dist;

    AKANTU_DEBUG_OUT();
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  Real value;
  BC::Axis incl_ax;
  Real center_coord;
  Real tang;
};

} // namespace akantu

#endif /* __AST_INCLINED_FLAT_DIRICHLET_HH__ */
