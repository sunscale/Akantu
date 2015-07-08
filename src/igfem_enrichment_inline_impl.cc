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
inline bool IGFEMEnrichment::isInside(const Vector<Real> & point, const ID & domain) {
  /// geometries[default_geometry];
  //  bool is_inside = true;
  Real tolerance = 1.e-14;
  Real radius = 0.401;

  /// @todo change to make this function generic!!
  Vector<Real> center(2);
  center(0) = 0.;
  center(1) = 0.;

  // if (std::abs(point(1)) >= point(0) + radius + tolerance)
  //   is_inside = false;

  return (point.distance(center) < radius + tolerance);
  ///return is_inside;
}

