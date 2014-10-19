/**
 * @file   data_accessor.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jun 16 2011
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  Interface of accessors for pack_unpack system
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


#ifndef __AKANTU_DATA_ACCESSOR_HH__
#define __AKANTU_DATA_ACCESSOR_HH__


/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "fe_engine.hh"
#include "communication_buffer.hh"
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

class DataAccessor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  DataAccessor();
  virtual ~DataAccessor();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /**
   * @brief get  the number of  data to exchange  for a given akantu::Element  and a
   * given akantu::SynchronizationTag
   */
  virtual UInt getNbDataForElements(__attribute__((unused)) const Array<Element> & elements,
				    __attribute__((unused)) SynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * @brief get  the number of  data to send  for a given
   * akantu::SynchronizationTag
   */
  virtual UInt getNbDataToPack(__attribute__((unused)) SynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * @brief get the number of data  to receive for a given
   * akantu::SynchronizationTag
   */
  virtual UInt getNbDataToUnpack(__attribute__((unused)) SynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * @brief   pack  the   data  for   a  given   akantu::Element  and   a  given
   * akantu::SynchronizationTag
   */
  virtual void packElementData(__attribute__((unused)) CommunicationBuffer & buffer,
			       __attribute__((unused)) const Array<Element> & element,
			       __attribute__((unused)) SynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * @brief   pack  the   data  for   a  given  index  and   a  given
   * akantu::SynchronizationTag
   */
  virtual void packData(__attribute__((unused)) CommunicationBuffer & buffer,
			__attribute__((unused)) const UInt index,
                        __attribute__((unused)) SynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * @brief   unpack  the   data  for   a  given   akantu::Element  and   a  given
   * akantu::SynchronizationTag
   */
  virtual void unpackElementData(__attribute__((unused)) CommunicationBuffer & buffer,
				 __attribute__((unused)) const Array<Element> & element,
				 __attribute__((unused)) SynchronizationTag tag) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /**
   * @brief   unpack  the   data  for   a  given  index  and   a  given
   * akantu::SynchronizationTag
   */
  virtual void unpackData(__attribute__((unused)) CommunicationBuffer & buffer,
                          __attribute__((unused)) const UInt index,
                          __attribute__((unused)) SynchronizationTag tag) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

public:
  template<typename T>
  static inline void packNodalDataHelper(const Array<T> & data,
					 CommunicationBuffer & buffer,
					 const Array<Element> & elements,
					 const Mesh & mesh) {
    packUnpackNodalDataHelper<T, true>(const_cast<Array<T> &>(data), 
				       buffer, 
				       elements, 
				       mesh);
  }

  template<typename T>
  static inline void unpackNodalDataHelper(Array<T> & data,
					   CommunicationBuffer & buffer,
					   const Array<Element> & elements,
					   const Mesh & mesh) {
    packUnpackNodalDataHelper<T, false>(data, 
					buffer, 
					elements, 
					mesh);
  }

  template<typename T, bool pack_helper>
  static inline void packUnpackNodalDataHelper(Array<T> & data,
					       CommunicationBuffer & buffer,
					       const Array<Element> & elements,
					       const Mesh & mesh);

  template<typename T, bool pack_helper>
  static inline void packUnpackElementalDataHelper(ElementTypeMapArray<T> & data_to_pack,
						   CommunicationBuffer & buffer,
						   const Array<Element> & element,
						   bool per_quadrature_point_data,
						   const FEEngine & fem);

  template<typename T>
  static inline void packElementalDataHelper(const ElementTypeMapArray<T> & data_to_pack,
					     CommunicationBuffer & buffer,
					     const Array<Element> & elements,
					     bool per_quadrature_point,
					     const FEEngine & fem) {
    packUnpackElementalDataHelper<T, true>(const_cast<ElementTypeMapArray<T> &>(data_to_pack),
					   buffer,
					   elements,
					   per_quadrature_point,
					   fem);
  }
  
  template<typename T>
  inline void unpackElementalDataHelper(ElementTypeMapArray<T> & data_to_unpack,
                                        CommunicationBuffer & buffer,
                                        const Array<Element> & elements,
                                        bool per_quadrature_point,
					const FEEngine & fem) {
    packUnpackElementalDataHelper<T, false>(data_to_unpack,
					    buffer,
					    elements,
					    per_quadrature_point,
					    fem);
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "data_accessor_inline_impl.cc"

// /// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const DataAccessor & _this)
// {
//   _this.printself(stream);
//   return stream;
// }

__END_AKANTU__

#endif /* __AKANTU_DATA_ACCESSOR_HH__ */
