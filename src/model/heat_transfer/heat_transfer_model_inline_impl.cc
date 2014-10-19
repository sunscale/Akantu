/**
 * @file   heat_transfer_model_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  Implementation of the inline functions of the HeatTransferModel class
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
inline FEEngine & HeatTransferModel::getFEEngineBoundary(std::string name) {
  return dynamic_cast<FEEngine &>(getFEEngineClassBoundary<MyFEEngineType>(name));
}

/* -------------------------------------------------------------------------- */
inline UInt HeatTransferModel::getNbDataToPack(SynchronizationTag tag) const{
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes = getFEEngine().getMesh().getNbNodes();

  switch(tag) {
  case _gst_htm_temperature:
  case _gst_htm_capacity: {
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
inline UInt HeatTransferModel::getNbDataToUnpack(SynchronizationTag tag) const{
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes = getFEEngine().getMesh().getNbNodes();

  switch(tag) {
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
					const UInt index,
					SynchronizationTag tag) const{
  AKANTU_DEBUG_IN();

  switch(tag) {
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

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::unpackData(CommunicationBuffer & buffer,
					  const UInt index,
					  SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  switch(tag) {
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

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline UInt HeatTransferModel::getNbDataForElements(const Array<Element> & elements,
                                                    SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = 0;
  Array<Element>::const_iterator<Element> it  = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

#ifndef AKANTU_NDEBUG
  size += elements.getSize() * spatial_dimension * sizeof(Real); /// position of the barycenter of the element (only for check)
  //  size += spatial_dimension * nb_nodes_per_element * sizeof(Real); /// position of the nodes of the element
#endif

  switch(tag) {
  case _gst_htm_capacity: {
    size += nb_nodes_per_element * sizeof(Real); // capacity vector
    break;
  }
  case _gst_htm_temperature: {
    size += nb_nodes_per_element * sizeof(Real); // temperature
    break;
  }
  case _gst_htm_gradient_temperature: {
    size += getNbQuadraturePoints(elements) * spatial_dimension * sizeof(Real); // temperature gradient
    size += nb_nodes_per_element * sizeof(Real); // nodal temperatures
    //    size += spatial_dimension * nb_nodes_per_element * sizeof(Real); // shape derivatives
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
inline void HeatTransferModel::packElementData(CommunicationBuffer & buffer,
                                               const Array<Element> & elements,
                                               SynchronizationTag tag) const {
#ifndef AKANTU_NDEBUG
  Array<Element>::const_iterator<Element> bit  = elements.begin();
  Array<Element>::const_iterator<Element> bend = elements.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;
    Vector<Real> barycenter(spatial_dimension);
    mesh.getBarycenter(element.element, element.type, barycenter.storage(), element.ghost_type);
    buffer << barycenter;
  }
  // packNodalDataHelper(mesh.getNodes(), buffer, elements);
#endif

  switch (tag){
  case _gst_htm_capacity: {
    packNodalDataHelper(*capacity_lumped, buffer, elements, mesh);
    break;
  }
  case _gst_htm_temperature: {
    packNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  case _gst_htm_gradient_temperature: {
    packElementalDataHelper(temperature_gradient, buffer, elements, true, getFEEngine());
    packNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::unpackElementData(CommunicationBuffer & buffer,
                                                 const Array<Element> & elements,
                                                 SynchronizationTag tag) {
#ifndef AKANTU_NDEBUG
  Array<Element>::const_iterator<Element> bit  = elements.begin();
  Array<Element>::const_iterator<Element> bend = elements.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;

    Vector<Real> barycenter_loc(spatial_dimension);
    mesh.getBarycenter(element.element, element.type, barycenter_loc.storage(), element.ghost_type);

    Vector<Real> barycenter(spatial_dimension);
    buffer >> barycenter;
    Real tolerance = 1e-15;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      if(!(std::abs(barycenter(i) - barycenter_loc(i)) <= tolerance))
	AKANTU_DEBUG_ERROR("Unpacking an unknown value for the element: "
			   << element
			   << "(barycenter[" << i << "] = " << barycenter_loc(i)
			   << " and buffer[" << i << "] = " << barycenter(i) << ") - tag: " << tag);
    }
  }

  // Vector<Real> coords(spatial_dimension);
  // Real * nodes = getFEEngine().getMesh().getNodes().storage();
  // for (UInt n = 0; n < nb_nodes_per_element; ++n) {
  //   buffer >> coords;
  //   UInt offset_conn = conn[el_offset + n];
  //   Real * coords_local = nodes+spatial_dimension*offset_conn;
  //   for (UInt i = 0; i < spatial_dimension; ++i) {
  //     if(!(std::abs(coords(i) - coords_local[i]) <= tolerance))
  //       AKANTU_EXCEPTION("Unpacking to wrong node for the element : "
  //       		 << element
  //       		 << "(coords[" << i << "] = " << coords_local[i]
  //       		 << " and buffer[" << i << "] = " << coords(i) << ")");
  //   }
  // }
#endif

  switch (tag){
  case _gst_htm_capacity: {
    unpackNodalDataHelper(*capacity_lumped, buffer, elements, mesh);
    break;
  }
  case _gst_htm_temperature: {
    unpackNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  case _gst_htm_gradient_temperature: {
    unpackElementalDataHelper(temperature_gradient, buffer, elements, true, getFEEngine());
    unpackNodalDataHelper(*temperature, buffer, elements, mesh);

    // //    Real tolerance = 1e-15;
    // if (!Math::are_vector_equal(spatial_dimension,gtemp.storage(),it_gtemp[element.element].storage())){
    //   Real dist = Math::distance_3d(gtemp.storage(), it_gtemp[element.element].storage());
    //   debug::debugger.getOutputStream().precision(20);
    //   std::stringstream temperatures_str;
    //   temperatures_str.precision(20);
    //   temperatures_str << std::scientific << "temperatures are ";
    //   for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    //     UInt offset_conn = conn[el_offset + n];
    //     temperatures_str << (*temperature)(offset_conn) << " ";
    //   }
    //   Array<Real>::matrix_iterator it_shaped =
    //     const_cast<Array<Real> &>(getFEEngine().getShapesDerivatives(element.type, ghost_type))
    //     .begin(nb_nodes_per_element,spatial_dimension);


    //   AKANTU_EXCEPTION("packed gradient do not match for element " << element.element << std::endl
    //     	       << "buffer is " << gtemp << " local is " << it_gtemp[element.element]
    //     	       << " dist is " << dist << std::endl
    //     	       << temperatures_str.str() << std::endl
    //     	       << std::scientific << std::setprecision(20)
    //     	       << " distant temperatures " << temp_nodes
    //     	       << "temperature gradient size " << temperature_gradient(element.type, ghost_type).getSize()
    //     	       << " number of ghost elements " << getFEEngine().getMesh().getNbElement(element.type,_ghost)
    //     	       << std::scientific << std::setprecision(20)
    //     	       << " shaped " << shaped
    //     	       << std::scientific << std::setprecision(20)
    //     	       << " local shaped " << it_shaped[element.element]);
    // }
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
