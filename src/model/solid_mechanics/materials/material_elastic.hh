/**
 * @file   material_elastic.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Nov 15 2015
 *
 * @brief  Material isotropic elastic
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "aka_common.hh"
#include "material_thermal.hh"
#include "plane_stress_toolbox.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_ELASTIC_HH__
#define __AKANTU_MATERIAL_ELASTIC_HH__

__BEGIN_AKANTU__

/**
 * Material elastic isotropic
 *
 * parameters in the material files :
 *   - E   : Young's modulus (default: 0)
 *   - nu  : Poisson's ratio (default: 1/2)
 *   - Plane_Stress : if 0: plane strain, else: plane stress (default: 0)
 */
template<UInt spatial_dimension>
class MaterialElastic : public PlaneStressToolbox< spatial_dimension,
                                                   MaterialThermal<spatial_dimension> > {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:
  typedef  PlaneStressToolbox< spatial_dimension,
                               MaterialThermal<spatial_dimension> > Parent;
public:

  MaterialElastic(SolidMechanicsModel & model, const ID & id = "");
  MaterialElastic(SolidMechanicsModel & model,
                  UInt dim,
                  const Mesh & mesh,
                  FEEngine & fe_engine,
                  const ID & id = "");

  virtual ~MaterialElastic() {}

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void initMaterial();

  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  virtual void computeTangentModuli(const ElementType & el_type,
				    Array<Real> & tangent_matrix,
				    GhostType ghost_type = _not_ghost);

  /// compute the elastic potential energy
  virtual void computePotentialEnergy(ElementType el_type,
				      GhostType ghost_type = _not_ghost);

  virtual void computePotentialEnergyByElement(ElementType type, UInt index,
					       Vector<Real> & epot_on_quad_points);

  /// compute the p-wave speed in the material
  virtual Real getPushWaveSpeed(const Element & element) const;

  /// compute the s-wave speed in the material
  virtual Real getShearWaveSpeed(const Element & element) const;

protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(const Matrix<Real> & grad_u,
				  Matrix<Real> & sigma,
				  const Real sigma_th = 0) const;

  /// compute the tangent stiffness matrix for an element
  inline void computeTangentModuliOnQuad(Matrix<Real> & tangent) const;

  /// recompute the lame coefficient if E or nu changes
  virtual void updateInternalParameters();

  static inline void computePotentialEnergyOnQuad(const Matrix<Real> & grad_u,
						  const Matrix<Real> & sigma,
						  Real & epot);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get first Lame constant
  AKANTU_GET_MACRO(Lambda, lambda, Real);

  /// get second Lame constant
  AKANTU_GET_MACRO(Mu, mu, Real);

  /// get bulk modulus
  AKANTU_GET_MACRO(Kappa, kpa, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// First Lamé coefficient
  Real lambda;

  /// Second Lamé coefficient (shear modulus)
  Real mu;

  /// Bulk modulus
  Real kpa;

};

__END_AKANTU__

#include "material_elastic_inline_impl.cc"

#endif /* __AKANTU_MATERIAL_ELASTIC_HH__ */
