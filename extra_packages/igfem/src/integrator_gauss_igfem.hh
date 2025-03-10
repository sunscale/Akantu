/**
 * @file   integrator_gauss_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Gauss integration facilities for IGFEM
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_INTEGRATOR_IGFEM_HH__
#define __AKANTU_INTEGRATOR_IGFEM_HH__

/* -------------------------------------------------------------------------- */
#include "integrator.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {

template <class IOF> class IntegratorGauss<_ek_igfem, IOF> : public Integrator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  IntegratorGauss(const Mesh & mesh, const ID & id = "integrator_gauss",
                  const MemoryID & memory_id = 0);

  virtual ~IntegratorGauss(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline void initIntegrator(const Array<Real> & nodes,
                             const ElementType & type,
                             const GhostType & ghost_type);

  /// precompute jacobians on elements of type "type"
  template <ElementType type>
  void precomputeJacobiansOnQuadraturePoints(const Array<Real> & nodes,
                                             const GhostType & ghost_type);

  /// integrate f on the element "elem" of type "type"
  template <ElementType type>
  inline void integrateOnElement(const Array<Real> & f, Real * intf,
                                 UInt nb_degree_of_freedom, const UInt elem,
                                 const GhostType & ghost_type) const;

  /// integrate f for all elements of type "type"
  template <ElementType type>
  void integrate(const Array<Real> & in_f, Array<Real> & intf,
                 UInt nb_degree_of_freedom, const GhostType & ghost_type,
                 const Array<UInt> & filter_elements) const;

  /// integrate one element scalar value on all elements of type "type"
  template <ElementType type>
  Real integrate(const Vector<Real> & in_f, UInt index,
                 const GhostType & ghost_type) const;

  /// integrate scalar field in_f
  template <ElementType type>
  Real integrate(const Array<Real> & in_f, const GhostType & ghost_type,
                 const Array<UInt> & filter_elements) const;

  /// integrate partially around a quadrature point (@f$ intf_q = f_q * J_q *
  /// w_q @f$)
  template <ElementType type>
  void integrateOnIntegrationPoints(const Array<Real> & in_f,
                                    Array<Real> & intf,
                                    UInt nb_degree_of_freedom,
                                    const GhostType & ghost_type,
                                    const Array<UInt> & filter_elements) const;
  /// return a vector with quadrature points natural coordinates
  template <ElementType type>
  const Matrix<Real> & getIntegrationPoints(const GhostType & ghost_type) const;

  /// return the number of quadrature points for a given element type
  template <ElementType type>
  inline UInt
  getNbIntegrationPoints(const GhostType & ghost_type = _not_ghost) const;

  /// compute the vector of quadrature points natural coordinates
  template <ElementType type>
  void computeQuadraturePoints(const GhostType & ghost_type);

  /// check that the jacobians are not negative
  template <ElementType type>
  void checkJacobians(const GhostType & ghost_type) const;

public:
  /// compute the jacobians on quad points for a given element
  template <ElementType type>
  void computeJacobianOnQuadPointsByElement(const Matrix<Real> & node_coords,
                                            Vector<Real> & jacobians);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  inline void integrate(Real * f, Real * jac, Real * inte,
                        UInt nb_degree_of_freedom,
                        UInt nb_quadrature_points) const;

private:
  ElementTypeMap<Matrix<Real>> quadrature_points;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "integrator_gauss_igfem_inline_impl.hh"

} // namespace akantu

#endif /*__AKANTU_INTEGRATOR_IGFEM_HH__ */
/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
