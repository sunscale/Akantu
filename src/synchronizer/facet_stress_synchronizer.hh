/**
 * @file   facet_stress_synchronizer.hh
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

#ifndef __AKANTU_FACET_STRESS_SYNCHRONIZER_HH__
#define __AKANTU_FACET_STRESS_SYNCHRONIZER_HH__

#include "distributed_synchronizer.hh"
#include "facet_synchronizer.hh"
#include "cohesive_element_inserter.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class FacetStressSynchronizer : public DistributedSynchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:

  FacetStressSynchronizer(Mesh & mesh_facets,
			  SynchronizerID id = "facet_stress_synchronizer",
			  MemoryID memory_id = 0);

  // virtual ~FacetStressSynchronizer();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  static FacetStressSynchronizer *
  createFacetStressSynchronizer(FacetSynchronizer & facet_synchronizer,
				Mesh & mesh_facets,
				SynchronizerID id = "facet_stress_synchronizer",
				MemoryID memory_id = 0);

  void updateFacetStressSynchronizer(const CohesiveElementInserter & inserter,
				     const ElementTypeMapArray<UInt> & rank_to_element,
				     DataAccessor & data_accessor);

protected:

  void updateElementList(Array<Element> * elements,
			 const CohesiveElementInserter & inserter,
			 const ElementTypeMapArray<UInt> & rank_to_element);

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

// #include "facet_stress_synchronizer_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_FACET_STRESS_SYNCHRONIZER_HH__ */
