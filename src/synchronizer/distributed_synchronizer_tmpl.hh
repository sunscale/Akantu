/**
 * @file   distributed_synchronizer_tmpl.hh
 *
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue May 07 2013
 * @date last modification: Mon Jun 09 2014
 *
 * @brief  Implementation of the templated function of the DistributedSynchronizer
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
#ifndef __AKANTU_DISTRIBUTED_SYNCHRONIZER_TMPL_HH__
#define __AKANTU_DISTRIBUTED_SYNCHRONIZER_TMPL_HH__

__BEGIN_AKANTU__

template<typename T>
void DistributedSynchronizer::fillTagBufferTemplated(const MeshData & mesh_data,
						     DynamicCommunicationBuffer * buffers,
						     const std::string & tag_name,
						     const ElementType & el_type,
						     const Array<UInt> & partition_num,
						     const CSR<UInt> & ghost_partition) {
  const Array<T> & data = mesh_data.getElementalDataArray<T>(tag_name, el_type);
  // Not possible to use the iterator because it potentially triggers the creation of complex
  // type templates (such as akantu::Vector< std::vector<Element> > which don't implement the right interface
  // (e.g. operator<< in that case).
  //typename Array<T>::template const_iterator< Vector<T> > data_it  = data.begin(data.getNbComponent());
  //typename Array<T>::template const_iterator< Vector<T> > data_end = data.end(data.getNbComponent());

  const T * data_it = data.storage();
  const T * data_end = data.storage() + data.getSize()*data.getNbComponent();
  const UInt * part = partition_num.storage();

  /// copying the data, element by element
  for (; data_it != data_end; ++part) {
    for(UInt j(0); j < data.getNbComponent(); ++j, ++data_it) {
      buffers[*part] << *data_it;
    }
  }

  data_it  = data.storage();
  /// copying the data for the ghost element
  for (UInt el(0); data_it != data_end; data_it+=data.getNbComponent(), ++el) {
    CSR<UInt>::const_iterator it = ghost_partition.begin(el);
    CSR<UInt>::const_iterator end = ghost_partition.end(el);
    for (;it != end; ++it) {
      UInt proc = *it;
      for(UInt j(0); j < data.getNbComponent(); ++j) {
	buffers[proc] << data_it[j];
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename BufferType>
void DistributedSynchronizer::populateMeshData(MeshData & mesh_data,
					       BufferType & buffer,
					       const std::string & tag_name,
					       const ElementType & el_type,
					       const MeshDataTypeCode & type_code,
					       UInt nb_component,
					       UInt nb_local_element,
					       UInt nb_ghost_element) {
#define AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA(r, extra_param, elem)	\
  case BOOST_PP_TUPLE_ELEM(2, 0, elem) : {				\
    populateMeshDataTemplated<BOOST_PP_TUPLE_ELEM(2, 1, elem)>(mesh_data, buffer, tag_name, el_type, nb_component, nb_local_element, nb_ghost_element); \
    break;								\
  }									\

  switch(type_code) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA, , AKANTU_MESH_DATA_TYPES)
  default : AKANTU_DEBUG_ERROR("Could not determine the type of tag" << tag_name << "!"); break;
  }
#undef AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA
}

/* -------------------------------------------------------------------------- */
template<typename T, typename BufferType>
void DistributedSynchronizer::populateMeshDataTemplated(MeshData & mesh_data,
							BufferType & buffer,
							const std::string & tag_name,
							const ElementType & el_type,
							UInt nb_component,
							UInt nb_local_element,
							UInt nb_ghost_element) {

  AKANTU_DEBUG_ASSERT(mesh.getNbElement(el_type) == nb_local_element,
		      "Did not got enought informations for the tag " << tag_name <<
		      " and the element type " << el_type << ":" << "_not_ghost." <<
		      " Got " << nb_local_element << " values, expected " << mesh.getNbElement(el_type));


  mesh_data.registerElementalData<T>(tag_name);
  Array<T> & data = mesh_data.getElementalDataArrayAlloc<T>(tag_name, el_type, _not_ghost, nb_component);
  data.resize(nb_local_element);
  /// unpacking the data, element by element
  for (UInt i(0); i < nb_local_element; ++i) {
    for(UInt j(0); j < nb_component; ++j) {
      buffer >> data(i,j);
    }
  }

  AKANTU_DEBUG_ASSERT(mesh.getNbElement(el_type, _ghost) == nb_ghost_element,
		      "Did not got enought informations for the tag " << tag_name <<
		      " and the element type " << el_type << ":" << "_ghost." <<
		      " Got " << nb_ghost_element << " values, expected " <<
		      mesh.getNbElement(el_type, _ghost));

  mesh_data.registerElementalData<T>(tag_name);
  Array<T> & data_ghost = mesh_data.getElementalDataArrayAlloc<T>(tag_name, el_type, _ghost, nb_component);
  data_ghost.resize(nb_ghost_element);

  /// unpacking the ghost data, element by element
  for (UInt j(0); j < nb_ghost_element; ++j) {
    for(UInt k(0); k < nb_component; ++k) {
      buffer >> data_ghost(j, k);
    }
  }
}

__END_AKANTU__


#endif /* __AKANTU_DISTRIBUTED_SYNCHRONIZER_TMPL_HH__ */
