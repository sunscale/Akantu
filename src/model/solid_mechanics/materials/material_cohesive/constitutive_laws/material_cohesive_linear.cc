/**
 * @file   material_cohesive_linear.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Thu Nov 19 2015
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
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
#include <algorithm>
#include <numeric>

/* -------------------------------------------------------------------------- */
#include "material_cohesive_linear.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "sparse_matrix.hh"
#include "dof_synchronizer.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialCohesiveLinear<spatial_dimension>::MaterialCohesiveLinear(SolidMechanicsModel & model,
                                                                  const ID & id) :
  MaterialCohesive(model,id),
  sigma_c_eff("sigma_c_eff", *this),
  delta_c_eff("delta_c_eff", *this),
  insertion_stress("insertion_stress", *this),
  opening_prec("opening_prec", *this),
  residual_sliding("residual_sliding", *this),
  friction_force("friction_force", *this),
  reduction_penalty("reduction_penalty", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("beta", beta, Real(0.),
                      _pat_parsable | _pat_readable,
                      "Beta parameter");

  this->registerParam("G_c", G_c, Real(0.),
                      _pat_parsable | _pat_readable,
                      "Mode I fracture energy");

  this->registerParam("penalty", penalty, Real(0.),
                      _pat_parsable | _pat_readable,
                      "Penalty coefficient");

  this->registerParam("volume_s", volume_s, Real(0.),
                      _pat_parsable | _pat_readable,
                      "Reference volume for sigma_c scaling");

  this->registerParam("m_s", m_s, Real(1.),
                      _pat_parsable | _pat_readable,
                      "Weibull exponent for sigma_c scaling");

  this->registerParam("kappa", kappa, Real(1.),
                      _pat_parsable | _pat_readable,
                      "Kappa parameter");

  this->registerParam("contact_after_breaking", contact_after_breaking, false,
		      _pat_parsable | _pat_readable,
		      "Activation of contact when the elements are fully damaged");

  this->registerParam("max_quad_stress_insertion", max_quad_stress_insertion, false,
		      _pat_parsable | _pat_readable,
		      "Insertion of cohesive element when stress is high enough just on one quadrature point");

  this->registerParam("friction", friction, false,
		      _pat_parsable | _pat_readable,
		      "Activation of friction in case of contact");

  this->registerParam("mu", mu_max, Real(0.),
                      _pat_parsable | _pat_readable,
                      "Maximum value of the friction coefficient");

  this->registerParam("penalty_for_friction", friction_penalty, Real(0.),
                      _pat_parsable | _pat_readable,
                      "Penalty parameter for the friction behavior");

  this->registerParam("recompute", recompute, false,
                      _pat_parsable | _pat_modifiable,
                      "recompute solution");

  use_previous_delta_max = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinear<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesive::initMaterial();

  /// compute scalars
  beta2_kappa2 = beta * beta/kappa/kappa;
  beta2_kappa  = beta * beta/kappa;

  if (Math::are_float_equal(beta, 0))
    beta2_inv = 0;
  else
    beta2_inv = 1./beta/beta;

  sigma_c_eff.initialize(1);
  delta_c_eff.initialize(1);
  insertion_stress.initialize(spatial_dimension);
  opening_prec.initialize(spatial_dimension);
  friction_force.initialize(spatial_dimension);
  residual_sliding.initialize(1);
  reduction_penalty.initialize(1);

  if (!Math::are_float_equal(delta_c, 0.))
    delta_c_eff.setDefaultValue(delta_c);
  else
    delta_c_eff.setDefaultValue(2 * G_c / sigma_c);

  if (model->getIsExtrinsic()) scaleInsertionTraction();
  residual_sliding.initializeHistory();
  opening_prec.initializeHistory();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinear<spatial_dimension>::scaleInsertionTraction() {
  AKANTU_DEBUG_IN();

  // do nothing if volume_s hasn't been specified by the user
  if (Math::are_float_equal(volume_s, 0.)) return;

  const Mesh & mesh_facets = model->getMeshFacets();
  const FEEngine & fe_engine = model->getFEEngine();
  const FEEngine & fe_engine_facet = model->getFEEngine("FacetsFEEngine");

  // loop over facet type
  Mesh::type_iterator first = mesh_facets.firstType(spatial_dimension - 1);
  Mesh::type_iterator last  = mesh_facets.lastType(spatial_dimension - 1);

  Real base_sigma_c = sigma_c;

  for(;first != last; ++first) {
    ElementType type_facet = *first;

    const Array< std::vector<Element> > & facet_to_element
      = mesh_facets.getElementToSubelement(type_facet);

    UInt nb_facet = facet_to_element.getSize();
    UInt nb_quad_per_facet = fe_engine_facet.getNbIntegrationPoints(type_facet);

    // iterator to modify sigma_c for all the quadrature points of a facet
    Array<Real>::vector_iterator sigma_c_iterator
      = sigma_c(type_facet).begin_reinterpret(nb_quad_per_facet, nb_facet);

    for (UInt f = 0; f < nb_facet; ++f, ++sigma_c_iterator) {

      const std::vector<Element> & element_list = facet_to_element(f);

      // compute bounding volume
      Real volume = 0;

      std::vector<Element>::const_iterator elem = element_list.begin();
      std::vector<Element>::const_iterator elem_end = element_list.end();

      for (; elem != elem_end; ++elem) {
        if (*elem == ElementNull) continue;

        // unit vector for integration in order to obtain the volume
        UInt nb_quadrature_points = fe_engine.getNbIntegrationPoints(elem->type);
        Vector<Real> unit_vector(nb_quadrature_points, 1);

        volume += fe_engine.integrate(unit_vector, elem->type,
                                      elem->element, elem->ghost_type);
      }

      // scale sigma_c
      *sigma_c_iterator -= base_sigma_c;
      *sigma_c_iterator *= std::pow(volume_s / volume, 1. / m_s);
      *sigma_c_iterator += base_sigma_c;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinear<spatial_dimension>::checkInsertion(bool check_only) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh_facets = model->getMeshFacets();
  CohesiveElementInserter & inserter = model->getElementInserter();

  Real tolerance = Math::getTolerance();

  Mesh::type_iterator it   = mesh_facets.firstType(spatial_dimension - 1);
  Mesh::type_iterator last = mesh_facets.lastType(spatial_dimension - 1);

  for (; it != last; ++it) {
    ElementType type_facet = *it;
    ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);
    const Array<bool> & facets_check = inserter.getCheckFacets(type_facet);
    Array<bool> & f_insertion = inserter.getInsertionFacets(type_facet);
    Array<UInt> & f_filter = facet_filter(type_facet);
    Array<Real> & sig_c_eff = sigma_c_eff(type_cohesive);
    Array<Real> & del_c = delta_c_eff(type_cohesive);
    Array<Real> & ins_stress = insertion_stress(type_cohesive);
    Array<Real> & trac_old = tractions_old(type_cohesive);
    Array<Real> & open_prec = opening_prec(type_cohesive);
    Array<Real> & res_sliding = residual_sliding(type_cohesive);
    Array<bool> & red_penalty = reduction_penalty(type_cohesive);
    const Array<Real> & f_stress = model->getStressOnFacets(type_facet);
    const Array<Real> & sigma_lim = sigma_c(type_facet);
    Real max_ratio = 0.;
    UInt index_f = 0;
    UInt index_filter = 0;
    UInt nn = 0;


    UInt nb_quad_facet = model->getFEEngine("FacetsFEEngine").getNbIntegrationPoints(type_facet);
    UInt nb_facet = f_filter.getSize();
    //  if (nb_facet == 0) continue;

    Array<Real>::const_iterator<Real> sigma_lim_it = sigma_lim.begin();

    Matrix<Real> stress_tmp(spatial_dimension, spatial_dimension);
    Matrix<Real> normal_traction(spatial_dimension, nb_quad_facet);
    Vector<Real> stress_check(nb_quad_facet);
    UInt sp2 = spatial_dimension * spatial_dimension;

    const Array<Real> & tangents = model->getTangents(type_facet);
    const Array<Real> & normals
      = model->getFEEngine("FacetsFEEngine").getNormalsOnIntegrationPoints(type_facet);
    Array<Real>::const_vector_iterator normal_begin = normals.begin(spatial_dimension);
    Array<Real>::const_vector_iterator tangent_begin = tangents.begin(tangents.getNbComponent());
    Array<Real>::const_matrix_iterator facet_stress_begin =
      f_stress.begin(spatial_dimension, spatial_dimension * 2);

    std::vector<Real> new_sigmas;
    std::vector< Vector<Real> > new_normal_traction;
    std::vector<Real> new_delta_c;

    // loop over each facet belonging to this material
    for (UInt f = 0; f < nb_facet; ++f, ++sigma_lim_it) {
      UInt facet = f_filter(f);
      // skip facets where check shouldn't be realized
      if (!facets_check(facet)) continue;

      // compute the effective norm on each quadrature point of the facet
      for (UInt q = 0; q < nb_quad_facet; ++q) {
	UInt current_quad = facet * nb_quad_facet + q;
	const Vector<Real> & normal = normal_begin[current_quad];
	const Vector<Real> & tangent = tangent_begin[current_quad];
	const Matrix<Real> & facet_stress_it = facet_stress_begin[current_quad];

	// compute average stress on the current quadrature point
	Matrix<Real> stress_1(facet_stress_it.storage(),
			      spatial_dimension,
			      spatial_dimension);

	Matrix<Real> stress_2(facet_stress_it.storage() + sp2,
			      spatial_dimension,
			      spatial_dimension);

	stress_tmp.copy(stress_1);
	stress_tmp += stress_2;
	stress_tmp /= 2.;

	Vector<Real> normal_traction_vec(normal_traction(q));

	// compute normal and effective stress
	stress_check(q) = computeEffectiveNorm(stress_tmp, normal, tangent,
					       normal_traction_vec);
      }

      // verify if the effective stress overcomes the threshold
      Real final_stress = stress_check.mean();
      if (max_quad_stress_insertion)
	final_stress = *std::max_element(stress_check.storage(),
					 stress_check.storage() + nb_quad_facet);

      if (final_stress > (*sigma_lim_it - tolerance)) {

        if (model->isExplicit()){
          f_insertion(facet) = true;

	  if (!check_only) {
	    // store the new cohesive material parameters for each quadrature point
	    for (UInt q = 0; q < nb_quad_facet; ++q) {
	      Real new_sigma = stress_check(q);
	      Vector<Real> normal_traction_vec(normal_traction(q));

	      if (spatial_dimension != 3)
		normal_traction_vec *= -1.;

	      new_sigmas.push_back(new_sigma);
	      new_normal_traction.push_back(normal_traction_vec);

	      Real new_delta;

	      // set delta_c in function of G_c or a given delta_c value
	      if (Math::are_float_equal(delta_c, 0.))
		new_delta = 2 * G_c / new_sigma;
	      else
		new_delta = (*sigma_lim_it) / new_sigma * delta_c;

	      new_delta_c.push_back(new_delta);
	    }
	  }

        }else{
          Real ratio = final_stress/(*sigma_lim_it);
          if (ratio > max_ratio){
            ++nn;
            max_ratio = ratio;
            index_f = f;
            index_filter = f_filter(f);
          }
        }
      }
    }

    /// Insertion of only 1 cohesive element in case of implicit approach. The one subjected to the highest stress.
    if (!model->isExplicit()){
      StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
      Array<Real> abs_max(comm.getNbProc());
      abs_max(comm.whoAmI()) = max_ratio;
      comm.allGather(abs_max.storage(), 1);

      Array<Real>::scalar_iterator it = std::max_element(abs_max.begin(), abs_max.end());
      Int pos = it - abs_max.begin();

      if (pos != comm.whoAmI()) {
      	AKANTU_DEBUG_OUT();
      	return;
      }

      if (nn) {
	f_insertion(index_filter) = true;

	if (!check_only) {
	  //  Array<Real>::iterator<Matrix<Real> > normal_traction_it =
	  //  normal_traction.begin_reinterpret(nb_quad_facet, spatial_dimension, nb_facet);
	  Array<Real>::const_iterator<Real> sigma_lim_it = sigma_lim.begin();

	  for (UInt q = 0; q < nb_quad_facet; ++q) {

	    //  Vector<Real> ins_s(normal_traction_it[index_f].storage() + q * spatial_dimension,
	    //            spatial_dimension);

	    Real new_sigma = (sigma_lim_it[index_f]);
            Vector<Real> normal_traction_vec(spatial_dimension, 0.0);

	    new_sigmas.push_back(new_sigma);
	    new_normal_traction.push_back(normal_traction_vec);

	    Real new_delta;

	    //set delta_c in function of G_c or a given delta_c value
	    if (!Math::are_float_equal(delta_c, 0.))
	      new_delta = delta_c;
	    else
	      new_delta = 2 * G_c / (new_sigma);

	    new_delta_c.push_back(new_delta);
	  }

	}
      }
    }

    // update material data for the new elements
    UInt old_nb_quad_points = sig_c_eff.getSize();
    UInt new_nb_quad_points = new_sigmas.size();
    sig_c_eff.resize(old_nb_quad_points + new_nb_quad_points);
    ins_stress.resize(old_nb_quad_points + new_nb_quad_points);
    trac_old.resize(old_nb_quad_points + new_nb_quad_points);
    del_c.resize(old_nb_quad_points + new_nb_quad_points);
    open_prec.resize(old_nb_quad_points + new_nb_quad_points);
    res_sliding.resize(old_nb_quad_points + new_nb_quad_points);
    red_penalty.resize(old_nb_quad_points + new_nb_quad_points);

    for (UInt q = 0; q < new_nb_quad_points; ++q) {
      sig_c_eff(old_nb_quad_points + q) = new_sigmas[q];
      del_c(old_nb_quad_points + q) = new_delta_c[q];
      res_sliding(old_nb_quad_points + q) = 0.;
      red_penalty(old_nb_quad_points + q) = false;
      for (UInt dim = 0; dim < spatial_dimension; ++dim) {
 	ins_stress(old_nb_quad_points + q, dim) = new_normal_traction[q](dim);
 	trac_old(old_nb_quad_points + q, dim) = new_normal_traction[q](dim);
 	open_prec(old_nb_quad_points + q, dim) = 0.;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinear<spatial_dimension>::computeTraction(const Array<Real> & normal,
                                                                ElementType el_type,
                                                                GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// define iterators
  Array<Real>::vector_iterator traction_it =
    tractions(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::vector_iterator opening_it =
    opening(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::vector_iterator opening_prec_it =
    opening_prec(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::vector_iterator opening_prec_prev_it =
    opening_prec.previous(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::vector_iterator contact_traction_it =
    contact_tractions(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::vector_iterator contact_opening_it =
    contact_opening(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::const_iterator< Vector<Real> > normal_it =
    normal.begin(spatial_dimension);

  Array<Real>::vector_iterator traction_end =
    tractions(el_type, ghost_type).end(spatial_dimension);

  Array<Real>::iterator<Real>sigma_c_it =
    sigma_c_eff(el_type, ghost_type).begin();

  Array<Real>::iterator<Real>delta_max_it =
    delta_max(el_type, ghost_type).begin();

  Array<Real>::iterator<Real>delta_max_prev_it =
    delta_max.previous(el_type, ghost_type).begin();

  Array<Real>::iterator<Real>delta_c_it =
    delta_c_eff(el_type, ghost_type).begin();

  Array<Real>::iterator<Real>damage_it =
    damage(el_type, ghost_type).begin();

  Array<Real>::vector_iterator insertion_stress_it =
    insertion_stress(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::iterator<Real>res_sliding_it =
    residual_sliding(el_type, ghost_type).begin();

  Array<Real>::iterator<Real>res_sliding_prev_it =
    residual_sliding.previous(el_type, ghost_type).begin();

  Array<Real>::vector_iterator friction_force_it =
    friction_force(el_type, ghost_type).begin(spatial_dimension);

  Array<bool>::iterator<bool> reduction_penalty_it =
    reduction_penalty(el_type, ghost_type).begin();

  Real * memory_space = new Real[2*spatial_dimension];
  Vector<Real> normal_opening(memory_space, spatial_dimension);
  Vector<Real> tangential_opening(memory_space + spatial_dimension,
                                  spatial_dimension);

  Vector<Real> open_prec(spatial_dimension);
  Vector<Real> open(spatial_dimension);
  Real max_delta_open = 0.0;
  Int nn = 0;

  /// loop on each quadrature point
  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++opening_prec_it, ++normal_it, ++sigma_c_it,
         ++delta_max_it, ++delta_c_it, ++damage_it, ++contact_traction_it,
         ++insertion_stress_it, ++contact_opening_it, ++delta_max_prev_it,
	 ++res_sliding_it, ++res_sliding_prev_it, ++opening_prec_prev_it,
	 ++friction_force_it, ++reduction_penalty_it) {

    if (!model->isExplicit())
      *delta_max_it = *delta_max_prev_it;

    // check if the local displacement increment is too large. If so,
    // reduce arbitrary it in order to prevent problems of
    // instability

    Real normal_opening_prec_norm = opening_prec_it->dot(*normal_it);

    open_prec = *opening_prec_it;
    Real open_prec_norm = open_prec.norm();
    open = *opening_it;
    Real open_norm = open.norm();
    Real delta_open = std::abs(open_norm - open_prec_norm);

    //    std::cout << "delta_open UPDATE RESIDUAL = " << delta_open << std::endl;
    max_delta_open = std::max(max_delta_open, delta_open);

    //    if (delta_open > 1.0e-01){
    //      *opening_it = *opening_prec_it + ((*opening_it - *opening_prec_it) / 50.0);
    //    }
    *opening_prec_it = *opening_it;

    Vector<Real> aa(spatial_dimension);
    aa = *opening_prec_it;
    //    std::cout << "opening_prec_it = " << aa(0) << ", " << aa(1) << std::endl;

    /// compute normal and tangential opening vectors
    Real normal_opening_norm = opening_it->dot(*normal_it);
    normal_opening  = (*normal_it);
    normal_opening *= normal_opening_norm;
    Vector<Real> perpend(spatial_dimension);
    perpend = *normal_it;
    //    std::cout << "normal = " << perpend(0) << perpend(1) << std::endl;
    Vector<Real> displ(spatial_dimension);
    displ = *opening_it;
    tangential_opening  = *opening_it;
    tangential_opening -= normal_opening;

    Real tangential_opening_norm = tangential_opening.norm();

    /**
     * compute effective opening displacement
     * @f$ \delta = \sqrt{
     * \frac{\beta^2}{\kappa^2} \Delta_t^2 + \Delta_n^2 } @f$
     */
    Real delta = tangential_opening_norm * tangential_opening_norm * beta2_kappa2;
    ++nn;
    //    std::cout << "normal opening of gp no. " << nn << " = " << normal_opening_norm << std::endl;

    bool penetration = normal_opening_norm < -Math::getTolerance();
    if (contact_after_breaking == false && Math::are_float_equal(*damage_it, 1.))
      penetration = false;

    /// if during the convergence loop a cohesive element continues to
    /// jumps from penetration to opening, and convergence is not
    /// reached, its penalty parameter will be reduced in the
    /// recomputation of the same incremental step
    if (!model->isExplicit() && !recompute)
      if ((normal_opening_prec_norm * normal_opening_norm) < -Math::getTolerance()){
	*reduction_penalty_it = true;
    
	//       std::cout << "recompute" << recompute << std::endl;
      }
    //    Vector<Real> friction_force(spatial_dimension);
    //    friction_force.clear();

    if (penetration) {
      Real current_penalty = 0.;
      if (recompute && *reduction_penalty_it){
	current_penalty = penalty / 1000.;
	//	std::cout << "recompute" << recompute << std::endl;
	//	std::cout << "penalty reduced" << std::endl;
      }
      else
	current_penalty = penalty;
      
      /// use penalty coefficient in case of penetration
      *contact_traction_it = normal_opening;
      *contact_traction_it *= current_penalty;
      *contact_opening_it = normal_opening;


      /* ----------------------------------------- */
      /// ADD FRICTION
      if (friction){
	Real damage = std::min(*delta_max_prev_it / *delta_c_it, Real(1.));
	Real mu = mu_max * damage;
	Vector<Real> check(spatial_dimension);
	check = *opening_prec_prev_it;
	Real normal_opening_prec_norm = opening_prec_prev_it->dot(*normal_it);
	Vector<Real> normal_opening_prec = (*normal_it);
	normal_opening_prec *= normal_opening_prec_norm;
	Real tau_max = mu * current_penalty * (normal_opening_prec.norm());
	Real delta_sliding_norm = std::abs(tangential_opening_norm - *res_sliding_prev_it);
	// tau is the norm of the friction force, acting tangentially to the surface
	Real tau = std::min(friction_penalty * delta_sliding_norm, tau_max);
	if ((tangential_opening_norm - *res_sliding_it) < -Math::getTolerance())
	  tau = -tau;
	// from tau get the x and y components of friction, to be added in the force vector
	Vector<Real> tangent(spatial_dimension);
	tangent = tangential_opening / tangential_opening_norm;
	*friction_force_it = tau * tangent;

	//	std::cout << "mu = " << mu << std::endl;
	//std::cout << "damage = " << damage << std::endl;
	//	std::cout << "opening_prev_prec_it = " << check(0) << ", " << check(1) << std::endl;
	//	std::cout << "normal opening prec norm = " << normal_opening_prec_norm << std::endl;
	//	std::cout << "normal opening prec = " << normal_opening_prec(0) << ", " << normal_opening_prec(1) << std::endl;
	//	std::cout << "tau_max = " << tau_max  << std::endl;

	//	std::cout << "friction force = " << friction_force(0) << ", " << friction_force(1) << std::endl;

	//update residual_sliding
	*res_sliding_it = tangential_opening_norm - (tau / friction_penalty);
      }
      /// end add friction
      /* ----------------------------------------- */

      /// don't consider penetration contribution for delta
      *opening_it = tangential_opening;
      normal_opening.clear();

    }
    else {
      delta += normal_opening_norm * normal_opening_norm;
      contact_traction_it->clear();
      contact_opening_it->clear();
      friction_force_it->clear();
    }

    delta = std::sqrt(delta);

    /// update maximum displacement and damage
    *delta_max_it = std::max(*delta_max_it, delta);
    *damage_it = std::min(*delta_max_it / *delta_c_it, Real(1.));

    /**
     * Compute traction @f$ \mathbf{T} = \left(
     * \frac{\beta^2}{\kappa} \Delta_t \mathbf{t} + \Delta_n
     * \mathbf{n} \right) \frac{\sigma_c}{\delta} \left( 1-
     * \frac{\delta}{\delta_c} \right)@f$
     */

    Vector<Real> forze(spatial_dimension);
    if (Math::are_float_equal(*damage_it, 1.))
      traction_it->clear();
    else if (Math::are_float_equal(*damage_it, 0.)) {
      if (penetration)
        traction_it->clear();
      else
        *traction_it = *insertion_stress_it;
    }
    else {
      *traction_it  = tangential_opening;
      *traction_it *= beta2_kappa;
      *traction_it += normal_opening;

      AKANTU_DEBUG_ASSERT(*delta_max_it != 0.,
                          "Division by zero, tolerance might be too low");

      *traction_it *= *sigma_c_it / *delta_max_it * (1. - *damage_it);

      forze = *traction_it;

      //      std::cout << "cohesive shear = " << forze(0) << ", " << forze(1) << std::endl;
    }

    if (friction)
      *traction_it += *friction_force_it;

    forze = *traction_it;

  }

  //  if (max_delta_open > 1.0e-01)
  //    std::cout << "max_delta opening = " << max_delta_open << std::endl;

  delete [] memory_space;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinear<spatial_dimension>::checkDeltaMax(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  
  /// This function set a predefined value to the parameter delta_max_prev of the
  /// elements that have been inserted in the last loading step for which convergence
  /// has not been reached. This is done before reducing the loading and re-doing the step.
  /// Otherwise, the updating of delta_max_prev would be done with reference to the
  /// non-convergent solution.

  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension,
					  ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension,
						ghost_type, _ek_cohesive);

  recompute = true;

  for(; it != last_type; ++it) {
    Array<UInt> & elem_filter = element_filter(*it, ghost_type);

    UInt nb_element = elem_filter.getSize();
    if (nb_element == 0) continue;

    ElementType el_type = *it;

    /// define iterators
    Array<Real>::iterator<Real>delta_max_it =
      delta_max(el_type, ghost_type).begin();

    Array<Real>::iterator<Real>delta_max_end =
      delta_max(el_type, ghost_type).end();

    Array<Real>::iterator<Real>delta_max_prev_it =
      delta_max.previous(el_type, ghost_type).begin();

    Array<Real>::iterator<Real>delta_c_it =
      delta_c_eff(el_type, ghost_type).begin();

    Array<Real>::vector_iterator opening_prec_it =
      opening_prec(el_type, ghost_type).begin(spatial_dimension);

    Array<Real>::vector_iterator opening_prec_prev_it =
      opening_prec.previous(el_type, ghost_type).begin(spatial_dimension);

    Array<Real>::iterator<Real>res_sliding_it =
      residual_sliding(el_type, ghost_type).begin();

    Array<Real>::iterator<Real>res_sliding_prev_it =
      residual_sliding.previous(el_type, ghost_type).begin();

    /// loop on each quadrature point
    for (; delta_max_it != delta_max_end;
         ++delta_max_it, ++delta_max_prev_it, ++delta_c_it,
	   ++opening_prec_it, ++opening_prec_prev_it,
	   ++res_sliding_it, ++res_sliding_prev_it) {

      if (*delta_max_prev_it == 0)
        /// elements inserted in the last incremental step, that did
        /// not converge
        *delta_max_it = *delta_c_it / 1000;
      else
        /// elements introduced in previous incremental steps, for
        /// which a correct value of delta_max_prev already exists
        *delta_max_it = *delta_max_prev_it;

      /// in case convergence is not reached, set opening_prec to the
      /// value referred to the previous incremental step
      *opening_prec_it = *opening_prec_prev_it;

      /// the same thing for the residual sliding for the friction law
      *res_sliding_it = *res_sliding_prev_it;
    }
  }
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinear<spatial_dimension>::resetVariables(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// This function set a predefined value to the parameter delta_max_prev of the
  /// elements that have been inserted in the last loading step for which convergence
  /// has not been reached. This is done before reducing the loading and re-doing the step.
  /// Otherwise, the updating of delta_max_prev would be done with reference to the
  /// non-convergent solution.

  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension,
					  ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension,
						ghost_type, _ek_cohesive);
  
  recompute = false;

  //  std::cout << "RESET VARIABLE" << std::endl;

  for(; it != last_type; ++it) {
    Array<UInt> & elem_filter = element_filter(*it, ghost_type);

    UInt nb_element = elem_filter.getSize();
    if (nb_element == 0) continue;

    ElementType el_type = *it;

    Array<bool>::iterator<bool> reduction_penalty_it =
      reduction_penalty(el_type, ghost_type).begin();

    Array<bool>::iterator<bool> reduction_penalty_end =
      reduction_penalty(el_type, ghost_type).end();

    /// loop on each quadrature point
    for (; reduction_penalty_it != reduction_penalty_end;
	 ++reduction_penalty_it) {
    
      *reduction_penalty_it = false;
    }
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveLinear<spatial_dimension>::computeTangentTraction(const ElementType & el_type,
                                                                       Array<Real> & tangent_matrix,
                                                                       const Array<Real> & normal,
                                                                       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// define iterators
  Array<Real>::matrix_iterator tangent_it =
    tangent_matrix.begin(spatial_dimension, spatial_dimension);

  Array<Real>::matrix_iterator tangent_end =
    tangent_matrix.end(spatial_dimension, spatial_dimension);

  Array<Real>::const_vector_iterator normal_it =
    normal.begin(spatial_dimension);

  Array<Real>::vector_iterator opening_it =
    opening(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::vector_iterator opening_prec_it =
    opening_prec.previous(el_type, ghost_type).begin(spatial_dimension);

  /// NB: delta_max_it points on delta_max_previous, i.e. the
  /// delta_max related to the solution of the previous incremental
  /// step
  Array<Real>::iterator<Real>delta_max_it =
    delta_max.previous(el_type, ghost_type).begin();

  Array<Real>::iterator<Real>sigma_c_it =
    sigma_c_eff(el_type, ghost_type).begin();

  Array<Real>::iterator<Real>delta_c_it =
    delta_c_eff(el_type, ghost_type).begin();

  Array<Real>::iterator<Real>damage_it =
    damage(el_type, ghost_type).begin();

  Array<Real>::iterator< Vector<Real> > contact_opening_it =
    contact_opening(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::iterator<Real>res_sliding_prev_it =
    residual_sliding.previous(el_type, ghost_type).begin();

  Array<bool>::iterator<bool> reduction_penalty_it =
    reduction_penalty(el_type, ghost_type).begin();

  Vector<Real> normal_opening(spatial_dimension);
  Vector<Real> tangential_opening(spatial_dimension);

  Vector<Real> open_prec(spatial_dimension);
  Vector<Real> open(spatial_dimension);
  Real max_delta_open = 0.0;

  Array<Real> tang_output(spatial_dimension,spatial_dimension);

  for (; tangent_it != tangent_end;
       ++tangent_it, ++normal_it, ++opening_it, ++opening_prec_it,
	 ++delta_max_it, ++sigma_c_it, ++delta_c_it, ++damage_it,
	 ++contact_opening_it, ++res_sliding_prev_it, ++reduction_penalty_it) {

    /// check if the local displacement increment is too large. In
    /// case, reduce arbitrary it in order to prevent problems of
    /// instability
    open_prec = *opening_prec_it;
    Real open_prec_norm = open_prec.norm();
    //    std::cout << "open_prec_norm = " << open_prec_norm << std::endl;

    /// During the update of the residual (called just before the
    /// assembly of the matrix) the interpenetrations are stored in
    /// the array "contact_opening", therefore, in the array "opening"
    /// there are only the tangential components.
    *opening_it += *contact_opening_it;

    open = *opening_it;
    Real open_norm = open.norm();
    //    std::cout << "open_norm = " << open_norm << std::endl;
    Real delta_open = std::abs(open_norm - open_prec_norm);

    max_delta_open = std::max(max_delta_open, delta_open);

    //    if (delta_open > 1.0e-01){
    //      *opening_it = *opening_prec_it + ((*opening_it - *opening_prec_it) / 50.0);
      //      std::cout << "delta_open TANGENT STIFFNESS = " << delta_open << std::endl;
    //    }

    /// compute normal and tangential opening vectors
    Real normal_opening_norm = opening_it->dot(*normal_it);
    normal_opening = (*normal_it);
    normal_opening *= normal_opening_norm;

    tangential_opening = *opening_it;
    tangential_opening -= normal_opening;

    Real tangential_opening_norm = tangential_opening.norm();
    bool penetration = normal_opening_norm < -Math::getTolerance();
    if (contact_after_breaking == false && Math::are_float_equal(*damage_it, 1.))
      penetration = false;

    Real derivative = 0; // derivative = d(t/delta)/ddelta
    Real t = 0;

    Real delta = tangential_opening_norm * tangential_opening_norm * beta2_kappa2;

    Matrix<Real> n_outer_n(spatial_dimension, spatial_dimension);
    n_outer_n.outerProduct(*normal_it, *normal_it);

    if (penetration){
      Real current_penalty = 0.;
      if (recompute && *reduction_penalty_it == 1)
	current_penalty = penalty / 1000.;
      else
	current_penalty = penalty;

      /// stiffness in compression given by the penalty parameter
      *tangent_it += n_outer_n;
      *tangent_it *= current_penalty;

      /* ------------------------------------- */
      /// ADD FRICTION
      if (friction){
	Real damage = std::min(*delta_max_it / *delta_c_it, Real(1.));
	Real mu = mu_max * damage;
	Real normal_opening_prec_norm = opening_prec_it->dot(*normal_it);
	Vector<Real> normal_opening_prec = (*normal_it);
	normal_opening_prec *= normal_opening_prec_norm;
	Real tau_max = mu * current_penalty * normal_opening_prec.norm();
	Real delta_sliding_norm = std::abs(tangential_opening_norm - *res_sliding_prev_it);
	// tau is the norm of the friction force, acting tangentially to the surface
	Real tau = std::min(friction_penalty * delta_sliding_norm, tau_max);

	if (tau < tau_max && tau_max > Math::getTolerance()){
	  Matrix<Real> I(spatial_dimension, spatial_dimension);
	  I.eye(1.);
	  Matrix<Real> nn(n_outer_n);
	  I -= nn;
	  *tangent_it += I * friction_penalty;
	}
      }
      /// end add friction
      /* ------------------------------------- */

      /// don't consider penetration contribution for delta
      //      std::cout << " AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA penetration !!!" << std::endl;

      *opening_it = tangential_opening;
      normal_opening_norm = opening_it->dot(*normal_it);
      normal_opening = (*normal_it);
      normal_opening *= normal_opening_norm;
      //      std::cout << "normal opening in case of penetration = " << normal_opening << std::endl;

    }
    else{
      delta += normal_opening_norm * normal_opening_norm;
    }

    delta = std::sqrt(delta);

    /// Delta has to be different from 0 to have finite values of tangential stiffness.
    /// At the element insertion, delta = 0. Therefore, a fictictious value is defined,
    /// for the evaluation of the first value of K.
    if (delta < Math::getTolerance())
      delta = (*delta_c_it)/1000.;

    //if (!penetration){
    if (delta >= *delta_max_it){
      if (delta <= *delta_c_it){
	derivative = -*sigma_c_it/(delta * delta);
	t = *sigma_c_it * (1 - delta / *delta_c_it);
      } else {
	//	std::cout << " BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB delta>delta_c !!!" << std::endl;
	derivative = 0.;
	t = 0.;
      }
    } else if (delta < *delta_max_it){
      Real tmax = *sigma_c_it * (1 - *delta_max_it / *delta_c_it);
      t = tmax / *delta_max_it * delta;
    }
      //    }


    /// computation of the derivative of the constitutive law (dT/ddelta)
    Matrix<Real> I(spatial_dimension, spatial_dimension);
    I.eye(beta2_kappa);
    Matrix<Real> nn(n_outer_n);
    nn *= (1 - beta2_kappa);
    nn += I;
    nn *= t/delta;
    Vector<Real> t_tilde(normal_opening);
    t_tilde *= (1 - beta2_kappa2);
    Vector<Real> mm(*opening_it);
    mm *= beta2_kappa2;
    t_tilde += mm;
    Vector<Real> t_hat(normal_opening);
    t_hat += beta2_kappa * tangential_opening;
    Matrix<Real> prov(spatial_dimension, spatial_dimension);
    prov.outerProduct(t_hat, t_tilde);
    prov *= derivative/delta;
    prov += nn;
    Matrix<Real> tmp(spatial_dimension, spatial_dimension);
    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = 0; j < spatial_dimension; ++j) {
	tmp(j,i) = prov(i,j);
      }
    }
    *tangent_it += tmp;

    /// check if the tangential stiffness matrix is symmetric
//    for (UInt h = 0; h < spatial_dimension; ++h){
//      for (UInt l = h; l < spatial_dimension; ++l){
//        if (l > h){
//          Real k_ls = (*tangent_it)[spatial_dimension*h+l];
//          Real k_us =  (*tangent_it)[spatial_dimension*l+h];
//          //          std::cout << "k_ls = " << k_ls << std::endl;
//          //          std::cout << "k_us = " << k_us << std::endl;
//          if (std::abs(k_ls) > 1e-13 && std::abs(k_us) > 1e-13){
//            Real error = std::abs((k_ls - k_us) / k_us);
//            if (error > 1e-10){
//	      std::cout << "non symmetric cohesive matrix" << std::endl;
//	      //  std::cout << "error " << error << std::endl;
//            }
//          }
//        }
//      }
//    }
  }

  //  if (max_delta_open > 1.0e-05)
  //    std::cout << "max_delta opening > 1e-05 in COMPUTE TANGENT STIFFNESS" << std::endl;
    
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */


INSTANTIATE_MATERIAL(MaterialCohesiveLinear);

__END_AKANTU__
