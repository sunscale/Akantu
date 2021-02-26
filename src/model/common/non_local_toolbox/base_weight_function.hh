/**
 * @file   base_weight_function.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Mon Aug 24 2015
 * @date last modification: Fri Dec 08 2017
 *
 * @brief  Base weight function for non local materials
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
#include "data_accessor.hh"
#include "model.hh"
#include "non_local_manager.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_BASE_WEIGHT_FUNCTION_HH_
#define AKANTU_BASE_WEIGHT_FUNCTION_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
/*  Normal weight function                                                    */
/* -------------------------------------------------------------------------- */
class BaseWeightFunction : public Parsable, public DataAccessor<Element> {
public:
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  BaseWeightFunction(NonLocalManager & manager,
                     const std::string & type = "base")
      : Parsable(ParserType::_weight_function, "weight_function:" + type),
        manager(manager), type(type),
        spatial_dimension(manager.getModel().getMesh().getSpatialDimension()) {
    this->registerParam("update_rate", update_rate, UInt(1), _pat_parsmod,
                        "Update frequency");
  }

  ~BaseWeightFunction() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  /// initialize the weight function
  virtual inline void init();

  /// update the internal parameters
  virtual void updateInternals(){};

  /* ------------------------------------------------------------------------ */
  /// set the non-local radius
  inline void setRadius(Real radius);

  /* ------------------------------------------------------------------------ */
  /// compute the weight for a given distance between two quadrature points
  inline Real operator()(Real r, const IntegrationPoint & q1,
                         const IntegrationPoint & q2) const;

  /// print function
  void printself(std::ostream & stream, int indent = 0) const override {
    std::string space;
    for (Int i = 0; i < indent; i++, space += AKANTU_INDENT) {
      ;
    }
    stream << space << "WeightFunction " << type << " [" << std::endl;
    Parsable::printself(stream, indent);
    stream << space << "]" << std::endl;
  }

  /* --------------------------------------------------------------------------
   */
  /* Accessors */
  /* --------------------------------------------------------------------------
   */

public:
  /// get the radius
  Real getRadius() const { return R; }
  /// get the update rate
  UInt getUpdateRate() const { return update_rate; }

public:
  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */

  UInt getNbData(const Array<Element> & /*elements*/,
                 const SynchronizationTag & /*tag*/) const override {
    return 0;
  }

  inline void packData(CommunicationBuffer & /*buffer*/,
                       const Array<Element> & /*element*/,
                       const SynchronizationTag & /*tag*/) const override {}

  inline void unpackData(CommunicationBuffer & /*buffer*/,
                         const Array<Element> & /*element*/,
                         const SynchronizationTag & /*tag*/) override {}

  /* ------------------------------------------------------------------------ */
  /* Accessors */
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

inline std::ostream & operator<<(std::ostream & stream,
                                 const BaseWeightFunction & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "base_weight_function_inline_impl.hh"

/* -------------------------------------------------------------------------- */
/* Include all other weight function types                                    */
/* -------------------------------------------------------------------------- */
#if defined(AKANTU_DAMAGE_NON_LOCAL)
#include "damaged_weight_function.hh"
#include "remove_damaged_weight_function.hh"
#include "remove_damaged_with_damage_rate_weight_function.hh"
#include "stress_based_weight_function.hh"
#endif
/* -------------------------------------------------------------------------- */

#endif /* AKANTU_BASE_WEIGHT_FUNCTION_HH_ */
