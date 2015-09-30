/**
 * @file   igfem_enrichment.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Implementation of IGFEM enrichment
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */
/* -------------------------------------------------------------------------- */
#include "igfem_enrichment.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
IGFEMEnrichment::IGFEMEnrichment(Mesh & mesh) : intersector_sphere(mesh) {
}

/* -------------------------------------------------------------------------- */
void IGFEMEnrichment::initialize(){
}

/* -------------------------------------------------------------------------- */
void IGFEMEnrichment::update(ID & domain) {
  if (domain == "") domain = default_geometry;
  Geometry & geometry = getGeometry(domain);
  intersector_sphere.buildIgfemMeshFromSpheres(geometry);
}
/* -------------------------------------------------------------------------- */
void IGFEMEnrichment::unRegisterGeometryObject(const ID & domain){
  
  GeometryMap::iterator it = geometries.find(domain);
  AKANTU_DEBUG_ASSERT(it != geometries.end(), "Geometry object with domain "
		      << domain << " was not found");
  geometries.erase(it);
  if (!geometries.empty())
    default_geometry = (*geometries.begin()).first;
}

/* -------------------------------------------------------------------------- */
void IGFEMEnrichment::registerGeometryObject(Geometry & geometry,
						    const ID & domain){
  if (geometries.size() == 0) default_geometry = domain;

#ifndef AKANTU_NDEBUG
  GeometryMap::iterator it = geometries.find(domain);
  AKANTU_DEBUG_ASSERT(it == geometries.end(), "Geometry object with domain "
		      << domain << " was already created");
#endif

  std::stringstream sstr; sstr << "geometry:" << domain;
  geometries[domain] = &geometry;
}

/* -------------------------------------------------------------------------- */
IGFEMEnrichment::Geometry & IGFEMEnrichment::getGeometry(ID & domain) const{
  AKANTU_DEBUG_IN();

  if (domain == "") domain = default_geometry;

  GeometryMap::const_iterator it = geometries.find(domain);
  AKANTU_DEBUG_ASSERT(it != geometries.end(), "The geometry " << domain << " is not registered");

  AKANTU_DEBUG_OUT();
  return *(it->second);
}

__END_AKANTU__
