/**
 * @file   igfem_enrichment_inline_impl.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  inline implementation of IGFEM enrichment
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

__END_AKANTU__

/// place here includes

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline bool IGFEMEnrichment::isInside(const Vector<Real> & point,
                                      ID domain) const {
  if (domain == "")
    domain = default_geometry;
  Geometry & spheres = this->getGeometry(domain);
  SK::Point_3 p(point(0), point(1), 0.);
  std::list<Spherical::Sphere_3>::const_iterator begin = spheres.begin();
  std::list<Spherical::Sphere_3>::const_iterator end = spheres.end();
  for (; begin != end; ++begin) {
    if (!(begin->has_on_unbounded_side(p)))
      return true;
  }
  return false;
}
