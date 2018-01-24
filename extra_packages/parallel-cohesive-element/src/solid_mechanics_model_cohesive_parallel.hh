/**
 * @file   solid_mechanics_model_cohesive_parallel.hh
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Include of the class members for parallelism
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------ */
/* Methods                                                                  */
/* ------------------------------------------------------------------------ */
public:
/// register the tags associated with the parallel synchronizer for
/// cohesive elements
void initParallel(MeshPartition * partition,
		  DataAccessor * data_accessor = NULL,
		  bool extrinsic = false);

protected:
void synchronizeGhostFacetsConnectivity();

void updateCohesiveSynchronizers();

/* ------------------------------------------------------------------------ */
/* Accessors                                                                */
/* ------------------------------------------------------------------------ */
public:

/// get cohesive elements synchronizer
AKANTU_GET_MACRO(CohesiveSynchronizer,
		 cohesive_distributed_synchronizer,
		 const DistributedSynchronizer *);

/* ------------------------------------------------------------------------ */
/* Data Accessor inherited members                                          */
/* ------------------------------------------------------------------------ */
public:

inline UInt getNbQuadsForFacetCheck(const Array<Element> & elements) const;

inline virtual UInt getNbDataForElements(const Array<Element> & elements,
					 SynchronizationTag tag) const;

inline virtual void packElementData(CommunicationBuffer & buffer,
				    const Array<Element> & elements,
				    SynchronizationTag tag) const;

inline virtual void unpackElementData(CommunicationBuffer & buffer,
				      const Array<Element> & elements,
				      SynchronizationTag tag);

template<typename T>
inline void packFacetStressDataHelper(const ElementTypeMapArray<T> & data_to_pack,
				      CommunicationBuffer & buffer,
				      const Array<Element> & elements) const;

template<typename T>
inline void unpackFacetStressDataHelper(ElementTypeMapArray<T> & data_to_unpack,
					CommunicationBuffer & buffer,
					const Array<Element> & elements) const;

template<typename T, bool pack_helper>
inline void packUnpackFacetStressDataHelper(ElementTypeMapArray<T> & data_to_pack,
					    CommunicationBuffer & buffer,
					    const Array<Element> & element) const;

/* ------------------------------------------------------------------------ */
/* Class Members                                                            */
/* ------------------------------------------------------------------------ */
private:
/// facet synchronizer
FacetSynchronizer * facet_synchronizer;

/// facet stress synchronizer
FacetStressSynchronizer * facet_stress_synchronizer;

/// cohesive elements synchronizer
DistributedSynchronizer * cohesive_distributed_synchronizer;

/// global connectivity
ElementTypeMapArray<UInt> * global_connectivity;
