/**
 * @file   base_weight_function.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Mon Aug 24 2015
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  Base weight function for non local materials
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "aka_common.hh"
#include "aka_types.hh"
#include "parsable.hh"
#include <cmath>
#if defined(AKANTU_DEBUG_TOOLS)
#include "aka_debug_tools.hh"
#include <string>
#endif
#include "non_local_manager.hh"

/* -------------------------------------------------------------------------- */
#include <vector>
#include "material_damage.hh"

#ifndef __AKANTU_BASE_WEIGHT_FUNCTION_HH__
#define __AKANTU_BASE_WEIGHT_FUNCTION_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
/*  Normal weight function                                                    */
/* -------------------------------------------------------------------------- */
class BaseWeightFunction : public Parsable {
public:
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  BaseWeightFunction(NonLocalManager & manager, const std::string & type = "base") :
    Parsable(_st_weight_function, "weight_function:" + type),
    manager(manager),
    type(type), 
    spatial_dimension(manager.getModel().getSpatialDimension()) {
    this->registerParam("update_rate"  , update_rate, 1U  ,
			_pat_parsmod, "Update frequency");    
  }

  virtual ~BaseWeightFunction() {}

  /* -------------------------------------------------------------------------- */
  /* Methods                                                                    */
  /* -------------------------------------------------------------------------- */
  /// initialize the weight function 
  virtual inline void init();

  /// update the internal parameters
  virtual void updateInternals() {};

  /* ------------------------------------------------------------------------ */
  /// set the non-local radius
  inline void setRadius(Real radius);

  /* ------------------------------------------------------------------------ */
  /// compute the weight for a given distance between two quadrature points
  inline Real operator()(Real r,
                         const __attribute__((unused)) IntegrationPoint & q1,
                         const __attribute__((unused)) IntegrationPoint & q2);

  /// print function
  void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
    stream << space << "WeightFunction " << type << " [" << std::endl;
    Parsable::printself(stream, indent);
    stream << space << "]" << std::endl;
  }

  /* -------------------------------------------------------------------------- */
  /* Accessors                                                                  */
  /* -------------------------------------------------------------------------- */

public:
  /// get the radius
  Real getRadius() { return R; }
  /// get the update rate
  UInt getUpdateRate() { return update_rate; }

public:

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */

  virtual UInt getNbDataForElements(__attribute__((unused)) const Array<Element> & elements,
                                    __attribute__((unused)) SynchronizationTag tag) const {
    return 0;
  }

  virtual inline void packElementData(__attribute__((unused)) CommunicationBuffer & buffer,
                                      __attribute__((unused)) const Array<Element> & elements,
                                      __attribute__((unused)) SynchronizationTag tag) const {}

  virtual inline void unpackElementData(__attribute__((unused)) CommunicationBuffer & buffer,
                                        __attribute__((unused)) const Array<Element> & elements,
                                        __attribute__((unused)) SynchronizationTag tag) {}

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Type, type, const ID &);

protected:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  /// reference to the non-local manager
  NonLocalManager & manager;

  /// the non-local radius
  Real R;

  /// the non-local radius squared
  Real R2;

  /// the update rate
  UInt update_rate;

  /// name of the type of weight function
  const std::string type;

  /// the spatial dimension 
  UInt spatial_dimension;
};

inline std::ostream & operator <<(std::ostream & stream,
                                  const BaseWeightFunction & _this)
{
  _this.printself(stream);
  return stream;
}


#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "base_weight_function_inline_impl.cc"
#endif

} // akantu
/* -------------------------------------------------------------------------- */
/* Include all other weight function types                                    */
/* -------------------------------------------------------------------------- */
#  include "damaged_weight_function.hh"
#  include "remove_damaged_weight_function.hh"
#  include "remove_damaged_with_damage_rate_weight_function.hh"
#  include "stress_based_weight_function.hh"



/* -------------------------------------------------------------------------- */


#endif /* __AKANTU_BASE_WEIGHT_FUNCTION_HH__ */
