/**
 * @file   dumper_field.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Tue Jan 19 2016
 *
 * @brief  Common interface for fields
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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

#ifndef __AKANTU_DUMPER_FIELD_HH__
#define __AKANTU_DUMPER_FIELD_HH__
/* -------------------------------------------------------------------------- */
#include "dumper_iohelper.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */
class FieldComputeProxy;
class FieldComputeBaseInterface;
class ComputeFunctorInterface;
class HomogenizerProxy;
/* -------------------------------------------------------------------------- */

/// Field interface
class Field {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  Field() = default;
  virtual ~Field() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
#ifdef AKANTU_USE_IOHELPER
  /// register this to the provided dumper
  virtual void registerToDumper(const std::string & id,
                                iohelper::Dumper & dumper) = 0;
#endif

  /// set the number of data per item (used for elements fields at the moment)
  virtual void setNbData(__attribute__((unused)) UInt nb_data) {
    AKANTU_TO_IMPLEMENT();
  };

  /// set the number of data per elem (used for elements fields at the moment)
  virtual void setNbDataPerElem(__attribute__((unused))
                                const ElementTypeMap<UInt> & nb_data) {
    AKANTU_TO_IMPLEMENT();
  };

  /// set the number of data per elem (used for elements fields at the moment)
  virtual void setNbDataPerElem(__attribute__((unused)) UInt nb_data) {
    AKANTU_TO_IMPLEMENT();
  };

  /// get the number of components of the hosted field
  virtual ElementTypeMap<UInt>
  getNbComponents(__attribute__((unused)) UInt dim = _all_dimensions,
                  __attribute__((unused)) GhostType ghost_type = _not_ghost,
                  __attribute__((unused)) ElementKind kind = _ek_not_defined) {
    throw;
  };

  /// for connection to a FieldCompute
  inline virtual Field * connect(__attribute__((unused))
                                 FieldComputeProxy & proxy) {
    throw;
  };

  /// for connection to a FieldCompute
  inline virtual ComputeFunctorInterface * connect(__attribute__((unused))
                                                   HomogenizerProxy & proxy) {
    throw;
  };

  /// check if the same quantity of data for all element types
  virtual void checkHomogeneity() = 0;

  /// return the dumper name
  std::string getGroupName() { return group_name; };

  /// return the id of the field
  std::string getID() { return field_id; };

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the flag to know if the field is homogeneous/contiguous
  virtual bool isHomogeneous() { return homogeneous; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the flag to know if it is homogeneous
  bool homogeneous{false};

  /// the name of the group it was associated to
  std::string group_name;

  /// the name of the dumper it was associated to
  std::string field_id;
};

/* -------------------------------------------------------------------------- */

__END_AKANTU_DUMPER__
} // akantu

#endif /* __AKANTU_DUMPER_FIELD_HH__ */
