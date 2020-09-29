/**
 * @file   igfem_enrichment.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  IGFEM enrichment: handles geometries of interfaces and advancement
 * with time
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_IGFEM_ENRICHMENT_HH_
#define AKANTU_IGFEM_ENRICHMENT_HH_

#include "mesh_igfem_spherical_growing_gel.hh"
#include "mesh_sphere_intersector.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

class IGFEMEnrichment {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  IGFEMEnrichment(Mesh & mesh);
  virtual ~IGFEMEnrichment(){};

private:
  typedef std::list<Spherical::Sphere_3> Geometry;
  typedef std::map<std::string, Geometry *> GeometryMap;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// create the mesh primitives
  void initialize();

  /// get a geometry from the geometry map
  virtual Geometry & getGeometry(ID & domain) const;

  /// detect the interface
  virtual void update(ID domain = "");

  /// remove geometry
  virtual void unRegisterGeometryObject(const ID & domain = "");

  /// insert new geometry
  virtual void registerGeometryObject(Geometry & geometry,
                                      const ID & domain = "");

  /// check if a point is in a given domain
  inline bool isInside(const Vector<Real> & point, ID domain = "") const;

  /// move the interface, in this case grow the gel pockets
  virtual void moveInterface(Real new_position, ID domain = "");

  /* --------------------------------------------------------------------------
   */
  /* Accessors */
  /* --------------------------------------------------------------------------
   */
public:
  UInt getNbStandardNodes() {
    return this->intersector_sphere.getNbStandardNodes();
  }

  UInt getNbEnrichedNodes() {
    return this->intersector_sphere.getNbEnrichedNodes();
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  MeshIgfemSphericalGrowingGel<2> intersector_sphere;

  GeometryMap geometries;

  /// default geometry object
  std::string default_geometry;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "igfem_enrichment_inline_impl.hh"

} // namespace akantu
/* -------------------------------------------------------------------------- */

#endif /* AKANTU_IGFEM_ENRICHMENT_HH_ */
