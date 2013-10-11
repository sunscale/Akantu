/**
 * @file   force_based_dirichlet.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Fri Oct 11 14:01:29 2013
 *
 * @brief  inclined dirichlet
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
#ifndef __AST_INCLINED_FLAT_DIRICHLET_HH__
#define __AST_INCLINED_FLAT_DIRICHLET_HH__

// akantu
#include "aka_common.hh"

// simtools
#include "ast_common.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
using namespace akantu;

/* -------------------------------------------------------------------------- */
class InclinedFlatDirichlet : public BC::Dirichlet::DirichletFunctor {
  
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  InclinedFlatDirichlet(Real val,
			BC::Axis ax,
			BC::Axis incl_ax,
			Real center_coord,
			Real tang) : 
    DirichletFunctor(ax),
    value(val),
    incl_ax(incl_ax),
    center_coord(center_coord),
    tang(tang) {};

  virtual ~InclinedFlatDirichlet() {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline void operator()(UInt node, 
			 Vector<bool> & flags, 
			 Vector<Real> & primal,
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

__END_SIMTOOLS__

#endif /* __AST_INCLINED_FLAT_DIRICHLET_HH__ */
