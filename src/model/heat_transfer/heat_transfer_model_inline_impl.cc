/**
 * @file   heat_transfer_model_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 *
 * @date creation: Fri Aug 20 2010
 * @date last modification: Tue Dec 08 2015
 *
 * @brief  Implementation of the inline functions of the HeatTransferModel class
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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

#ifndef __AKANTU_HEAT_TRANSFER_MODEL_INLINE_IMPL_CC__
#define __AKANTU_HEAT_TRANSFER_MODEL_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt HeatTransferModel::getNbData(const Array<UInt> & indexes,
                                         const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes = indexes.size();

  switch (tag) {
  case _gst_htm_capacity:
  case _gst_htm_temperature: {
    size += nb_nodes * sizeof(Real);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::packData(CommunicationBuffer & buffer,
                                        const Array<UInt> & indexes,
                                        const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  for (auto index : indexes) {
    switch (tag) {
    case _gst_htm_capacity:
      buffer << (*capacity_lumped)(index);
      break;
    case _gst_htm_temperature: {
      buffer << (*temperature)(index);
      break;
    }
    default: {
      AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
    }
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::unpackData(CommunicationBuffer & buffer,
                                          const Array<UInt> & indexes,
                                          const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  for (auto index : indexes) {
    switch (tag) {
    case _gst_htm_capacity: {
      buffer >> (*capacity_lumped)(index);
      break;
    }
    case _gst_htm_temperature: {
      buffer >> (*temperature)(index);
      break;
    }
    default: {
      AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
    }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline UInt HeatTransferModel::getNbData(const Array<Element> & elements,
                                         const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = 0;
  Array<Element>::const_iterator<Element> it = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

  switch (tag) {
  case _gst_htm_capacity: {
    size += nb_nodes_per_element * sizeof(Real); // capacity vector
    break;
  }
  case _gst_htm_temperature: {
    size += nb_nodes_per_element * sizeof(Real); // temperature
    break;
  }
  case _gst_htm_gradient_temperature: {
    // temperature gradient
    size += getNbIntegrationPoints(elements) * spatial_dimension * sizeof(Real);
    size += nb_nodes_per_element * sizeof(Real); // nodal temperatures
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::packData(CommunicationBuffer & buffer,
                                        const Array<Element> & elements,
                                        const SynchronizationTag & tag) const {
  switch (tag) {
  case _gst_htm_capacity: {
    packNodalDataHelper(*capacity_lumped, buffer, elements, mesh);
    break;
  }
  case _gst_htm_temperature: {
    packNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  case _gst_htm_gradient_temperature: {
    packElementalDataHelper(temperature_gradient, buffer, elements, true,
                            getFEEngine());
    packNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::unpackData(CommunicationBuffer & buffer,
                                          const Array<Element> & elements,
                                          const SynchronizationTag & tag) {
  switch (tag) {
  case _gst_htm_capacity: {
    unpackNodalDataHelper(*capacity_lumped, buffer, elements, mesh);
    break;
  }
  case _gst_htm_temperature: {
    unpackNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  case _gst_htm_gradient_temperature: {
    unpackElementalDataHelper(temperature_gradient, buffer, elements, true,
                              getFEEngine());
    unpackNodalDataHelper(*temperature, buffer, elements, mesh);

    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* __AKANTU_HEAT_TRANSFER_MODEL_INLINE_IMPL_CC__ */
