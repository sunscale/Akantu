/**
 * @file   facet_stress_synchronizer.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 *
 * @brief  Stress check on facets synchronizer
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

#include "facet_stress_synchronizer.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
FacetStressSynchronizer::FacetStressSynchronizer(Mesh & mesh_facets,
						 SynchronizerID id,
						 MemoryID memory_id) :
  DistributedSynchronizer(mesh_facets, id, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
FacetStressSynchronizer * FacetStressSynchronizer::
createFacetStressSynchronizer(FacetSynchronizer & facet_synchronizer,
			      Mesh & mesh_facets,
			      SynchronizerID id,
			      MemoryID memory_id) {
  AKANTU_DEBUG_IN();

  FacetStressSynchronizer & fs_synchronizer =
    *(new FacetStressSynchronizer(mesh_facets, id, memory_id));

  for (UInt p = 0; p < fs_synchronizer.nb_proc; ++p) {
    fs_synchronizer.send_element[p].copy(facet_synchronizer.recv_element[p]);
    fs_synchronizer.recv_element[p].copy(facet_synchronizer.send_element[p]);
  }

  AKANTU_DEBUG_OUT();
  return &fs_synchronizer;
}

/* -------------------------------------------------------------------------- */
void FacetStressSynchronizer::updateFacetStressSynchronizer(const CohesiveElementInserter & inserter,
							    const ElementTypeMapArray<UInt> & rank_to_element,
							    DataAccessor & data_accessor) {
  AKANTU_DEBUG_IN();

  updateElementList(send_element, inserter, rank_to_element);
  updateElementList(recv_element, inserter, rank_to_element);

  std::map<SynchronizationTag, Communication>::iterator it = communications.begin();
  std::map<SynchronizationTag, Communication>::iterator end = communications.end();

  for (; it != end; ++it) {
    SynchronizationTag tag = it->first;
    computeBufferSize(data_accessor, tag);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FacetStressSynchronizer::updateElementList(Array<Element> * elements,
						const CohesiveElementInserter & inserter,
						const ElementTypeMapArray<UInt> & rank_to_element) {
  AKANTU_DEBUG_IN();

  for (UInt p = 0; p < nb_proc; ++p) {

    ElementType current_element_type = _not_defined;
    GhostType current_ghost_type = _casper;
    const Array<bool> * f_check = NULL;
    const Array< std::vector<Element> > * element_to_facet = NULL;

    UInt nb_element = 0;

    Array<Element>::iterator<Element> it_new  = elements[p].begin();
    Array<Element>::iterator<Element> it  = elements[p].begin();
    Array<Element>::iterator<Element> end = elements[p].end();
    for (; it != end; ++it) {
      const Element & el = *it;

      if(el.type != current_element_type || el.ghost_type != current_ghost_type) {
	current_element_type = el.type;
	current_ghost_type   = el.ghost_type;

	element_to_facet = & mesh.getElementToSubelement(el.type, el.ghost_type);

	f_check = & inserter.getCheckFacets(el.type, el.ghost_type);
      }

      if ( (*f_check)(el.element) ) {

	const Element & el_0 = (*element_to_facet)(el.element)[0];
	const Array<UInt> & rank_el_0 = rank_to_element(el_0.type, el_0.ghost_type);

	const Element & el_1 = (*element_to_facet)(el.element)[1];
	const Array<UInt> & rank_el_1 = rank_to_element(el_1.type, el_1.ghost_type);

	if ( (rank_el_0(el_0.element) == rank && rank_el_1(el_1.element) == p) ||
	     (rank_el_0(el_0.element) == p && rank_el_1(el_1.element) == rank) ) {
	  *it_new = el;
	  ++it_new;
	  ++nb_element;
	}
      }
    }
    elements[p].resize(nb_element);
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
