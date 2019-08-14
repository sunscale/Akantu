/**
 * @file   random_internal_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Random internal material parameter
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_random_generator.hh"
#include "internal_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_RANDOM_INTERNAL_FIELD_HH__
#define __AKANTU_RANDOM_INTERNAL_FIELD_HH__

namespace akantu {

/**
 * class for the internal fields of materials with a random
 * distribution
 */
template <typename T, template <typename> class BaseField = InternalField,
          template <typename> class Generator = RandomGenerator>
class RandomInternalField : public BaseField<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  RandomInternalField(const ID & id, Material & material);

  ~RandomInternalField() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  RandomInternalField operator=(const RandomInternalField &) = delete;

public:
  AKANTU_GET_MACRO(RandomParameter, random_parameter, const RandomParameter<T>);

  /// initialize the field to a given number of component
  void initialize(UInt nb_component) override;

  /// set the field to a given value
  void setDefaultValue(const T & value) override;

  /// set the specified random distribution to a given parameter
  void setRandomDistribution(const RandomParameter<T> & param);

  /// print the content
  void printself(std::ostream & stream, int indent = 0) const override;

protected:
  void setArrayValues(T * begin, T * end) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  inline operator Real() const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// random parameter containing the distribution and base value
  RandomParameter<T> random_parameter;
};

/// standard output stream operator
template <typename T>
inline std::ostream & operator<<(std::ostream & stream,
                                 const RandomInternalField<T> & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* __AKANTU_RANDOM_INTERNAL_FIELD_HH__ */
