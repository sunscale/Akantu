/**
 * @file   base_weight_function.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Fri Apr 13 2012
 * @date last modification: Thu Jun 05 2014
 *
 * @brief Base weight function for non local materials
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
#include "aka_common.hh"
#include "aka_types.hh"
#include "solid_mechanics_model.hh"
#include "parsable.hh"
#include <cmath>
#if defined(AKANTU_DEBUG_TOOLS)
#include "aka_debug_tools.hh"
#include <string>
#endif

#include "material_list.hh"

/* -------------------------------------------------------------------------- */
#include <vector>
#include "material_damage.hh"

#ifndef __AKANTU_BASE_WEIGHT_FUNCTION_HH__
#define __AKANTU_BASE_WEIGHT_FUNCTION_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/*  Normal weight function                                                    */
/* -------------------------------------------------------------------------- */
class BaseWeightFunction : public Parsable {
public:
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  BaseWeightFunction(Material & material, const std::string & type = "base") :
    Parsable(_st_non_local, "weight_function:" + type), material(material), type(type), spatial_dimension(material.getSpatialDimension()) {
    this->registerParam("radius"       , R             , 100.,
			_pat_parsable | _pat_readable  , "Non local radius");
    this->registerParam("update_rate"  , update_rate, 0U  ,
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
  /// function for optimization to preselect types for the
  /// ElementTypeMapArrays (this function is currently not used)
  inline void selectType(__attribute__((unused)) ElementType type1,
                         __attribute__((unused)) GhostType   ghost_type1,
                         __attribute__((unused)) ElementType type2,
                         __attribute__((unused)) GhostType   ghost_type2) {
  }

  /* ------------------------------------------------------------------------ */
  /// compute the weight for a given distance between two quadrature points
  inline Real operator()(Real r,
                         const __attribute__((unused)) QuadraturePoint & q1,
                         const __attribute__((unused)) QuadraturePoint & q2);

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
  /// the material for which the weights are computed
  Material & material;

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

__END_AKANTU__

/* -------------------------------------------------------------------------- */


#endif /* __AKANTU_BASE_WEIGHT_FUNCTION_HH__ */
