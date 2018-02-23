/**
 * @file   ntn_base_contact_inline_impl.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  ntn base contact inline functions
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
inline UInt
NTNBaseContact::getNbDataForElements(const Array<Element> & elements,
                                     SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt spatial_dimension = this->model.getSpatialDimension();

  UInt nb_nodes = 0;
  Array<Element>::const_iterator<Element> it = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    nb_nodes += Mesh::getNbNodesPerElement(el.type);
  }

  switch (tag) {
  case _gst_cf_nodal: {
    size += nb_nodes * spatial_dimension * sizeof(Real) *
            3; // disp, vel and cur_pos
    break;
  }
  case _gst_cf_incr: {
    size += nb_nodes * spatial_dimension * sizeof(Real) * 1;
    break;
  }
  default: {}
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void NTNBaseContact::packElementData(CommunicationBuffer & buffer,
                                            const Array<Element> & elements,
                                            SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  switch (tag) {

  case _gst_cf_nodal: {
    DataAccessor::packNodalDataHelper(this->model.getDisplacement(), buffer,
                                      elements, this->model.getMesh());
    DataAccessor::packNodalDataHelper(this->model.getCurrentPosition(), buffer,
                                      elements, this->model.getMesh());
    DataAccessor::packNodalDataHelper(this->model.getVelocity(), buffer,
                                      elements, this->model.getMesh());
    break;
  }
  case _gst_cf_incr: {
    DataAccessor::packNodalDataHelper(this->model.getIncrement(), buffer,
                                      elements, this->model.getMesh());
    break;
  }
  default: {}
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void NTNBaseContact::unpackElementData(CommunicationBuffer & buffer,
                                              const Array<Element> & elements,
                                              SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  switch (tag) {

  case _gst_cf_nodal: {
    DataAccessor::unpackNodalDataHelper(this->model.getDisplacement(), buffer,
                                        elements, this->model.getMesh());
    DataAccessor::unpackNodalDataHelper(
        const_cast<Array<Real> &>(this->model.getCurrentPosition()), buffer,
        elements, this->model.getMesh());
    DataAccessor::unpackNodalDataHelper(this->model.getVelocity(), buffer,
                                        elements, this->model.getMesh());
    break;
  }
  case _gst_cf_incr: {
    DataAccessor::unpackNodalDataHelper(this->model.getIncrement(), buffer,
                                        elements, this->model.getMesh());
    break;
  }
  default: {}
  }

  AKANTU_DEBUG_OUT();
}
