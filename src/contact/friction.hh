/**
 * @file   friction.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Fri Jan 04 2013
 *
 * @brief  contact friction classes
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

#ifndef __AKANTU_FRICTION_HH__
#define __AKANTU_FRICTION_HH__

#include <iostream>
#include <cmath>


#include "aka_common.hh"

__BEGIN_AKANTU__

using std::cout;
using std::endl;


enum Friction_type { Coulomb_t, Prakash_Clifton_t, Rate_state_t, Rice_t };


// base class for all contract friction concrete classes,
// used so that friction classes can be stored in a container
// if needed
struct Contact_friction_base {};

template <int>
class Contact_friction;

// partial template specialization for Coulomb friction
template <>
class Contact_friction<Coulomb_t> : public Contact_friction_base {
    
    typedef Contact_friction_base base_type;
    typedef Real value_type;
    
protected:
    value_type mu_;  // coefficient of friction
    
public:
    
    Contact_friction() : base_type(), mu_() {}
    Contact_friction(value_type mu) : base_type(), mu_(mu) {}
    
    void setFriction(value_type mu)
    { mu_ = mu; }
    
    value_type friction() const
    { return mu_; }
    
    value_type computeForceMagnitude(value_type Pn)
    { return mu_*Pn; }

    template <class vector_type>
    value_type computeForceMagnitude(const vector_type& P) {
        
        value_type norm = value_type();
        for (int i=0; i<P.size(); ++i)
            norm += pow(mu_*P[i], 2.);
        return sqrt(norm);
    }
};




__END_AKANTU__

#endif /* __AKANTU_FRICTION_HH__ */
