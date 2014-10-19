/**
 * @file   resolution_augmented_lagrangian.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Mon Sep 15 2014
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  contact resolution classes
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_RESOLUTION_AUGMENTED_LAGRANGIAN_HH__
#define __AKANTU_RESOLUTION_AUGMENTED_LAGRANGIAN_HH__

#include <map>

#include "search.hh"
#include "aka_common.hh"
#include "contact_common.hh"
#include "aka_array.hh"
#include "dumpable_inline_impl.hh"
//#include "dumper_iohelper_tmpl.hh"
#include "solid_mechanics_model_element.hh"
#include "aka_optimize.hh"
#include "parsable.hh"

#define __LOC__                                                                \
  std::cout << "File " << __FILE__ << ", line " << __LINE__ << std::endl;

__BEGIN_AKANTU__

template <int Dim>
class ContactResolution<Dim, _static,
                        _augmented_lagrangian> : public ContactParameters,
                                                 public Parsable {
  typedef SolidMechanicsModel model_type;
  typedef ModelElement<model_type> element_type;

  typedef std::map<UInt, element_type> slave_master_map;
  typedef std::map<UInt, Real> real_map;
  typedef typename real_map::iterator real_iterator;
  typedef std::map<UInt, Real> gap_map;

  typedef Point<Dim> point_type;
  typedef array::Array<1, Real> vector_type;
  typedef array::Array<2, Real> matrix_type;

  model_type &model_;
  slave_master_map sm_;
  real_map areas_, gaps_;
  real_map multipliers_, penalty_;
  size_t uiter_, niter_;
  Array<Real> multiplier_dumper_, pressure_dumper_;
  std::map<UInt, bool> contact_status_;
  std::map<UInt, int> status_change_;

protected:
  //! Initialize function that's called after the derived Contact's constructor
  virtual void initialize();

public:
  //! Parameter constructor
  ContactResolution(model_type &m);

  //! Add slave
  void addSlave(UInt s) { sm_[s] = element_type(); }

  //! Add slave-master pair
  void addPair(UInt s, element_type el) { sm_[s] = el; }

  //! Add area to a slave node
  void addArea(UInt n, Real a) {
    if (a != 0.)
      areas_[n] = a;
  }

  //! Implementation of the contact step
  template <ContactImplementationMethod i>
  void solveContactStep(SearchBase *sp) {
    solveContactStepImpl(sp, Int2Type<i>());
  }

  //! Get contact force (sum of Lagrange multipliers)
  Real getForce() {
    Real f = 0.;
    for (auto v : multipliers_)
      f += v.second;
    return f;
  }

  //! Get maximum value of contact pressure
  Real getMaxPressure() {
    Real p = 0.;
    for (auto v : multipliers_) {
      auto slave = v.first;
      auto it = areas_.find(slave);
      if (it != areas_.end())
        p = std::max(p, v.second / it->second); // max(p, lambda/area)
    }
    return p;
  }

private:
  //! Compute penalty parameter values for each slave node automatically
  void getPenaltyValues();

  //! Add stiffness matrix and force vector contributions due to contact
  bool computeTangentAndResidual(real_map &lambda_new, SearchBase *cp,
                                 Real &error, Int2Type<_uzawa>);
  bool computeTangentAndResidual(Array<Real> &solution, Array<Real> &rhs,
                                 SearchBase *cp, Real &error,
                                 Int2Type<_generalized_newton>);

  //! Implementation of the contact step using Uzawa's method
  void solveContactStepImpl(SearchBase *sp, Int2Type<_uzawa>);

  //! Implementation of the contact step using Uzawa's method
  void solveContactStepImpl(SearchBase *sp, Int2Type<_generalized_newton>);

protected:
  //! Dump of Paraview files including contact information
  void dump();

  //! Provide standard output of resolution object
  template <int D>
  friend std::ostream &operator<<(
      std::ostream &os,
      const ContactResolution<D, _static, _augmented_lagrangian> &);
};

__END_AKANTU__

#endif /* __AKANTU_RESOLUTION_AUGMENTED_LAGRANGIAN_HH__ */
