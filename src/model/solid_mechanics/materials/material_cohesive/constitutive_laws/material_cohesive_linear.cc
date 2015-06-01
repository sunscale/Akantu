/**
 * @file   material_cohesive_linear.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Thu Aug 07 2014
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
  insertion_stress("insertion_stress", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("beta"   , beta   , 0. ,
                      _pat_parsable | _pat_readable,
                      "Beta parameter"         );

  this->registerParam("G_c"   , G_c   , 0. ,
                      _pat_parsable | _pat_readable,
                      "Mode I fracture energy" );

  this->registerParam("penalty", penalty, 0. ,
                      _pat_parsable | _pat_readable,
                      "Penalty coefficient"    );

  this->registerParam("volume_s", volume_s, 0. ,
                      _pat_parsable | _pat_readable,
                      "Reference volume for sigma_c scaling");

  this->registerParam("m_s", m_s, 1. ,
                      _pat_parsable | _pat_readable,
                      "Weibull exponent for sigma_c scaling");

  this->registerParam("kappa"  , kappa  , 1. ,
                      _pat_parsable | _pat_readable,
                      "Kappa parameter");

  //  if (!model->isExplicit())
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

  if (!Math::are_float_equal(delta_c, 0.))
    delta_c_eff.setDefaultValue(delta_c);
  else
    delta_c_eff.setDefaultValue(2 * G_c / sigma_c);

  if (model->getIsExtrinsic()) scaleInsertionTraction();

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
    UInt nb_quad_per_facet = fe_engine_facet.getNbQuadraturePoints(type_facet);

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
        UInt nb_quadrature_points = fe_engine.getNbQuadraturePoints(elem->type);
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
    const Array<Real> & f_stress = model->getStressOnFacets(type_facet);
    const Array<Real> & sigma_lim = sigma_c(type_facet);
    Real max_ratio = 0.;
    UInt index_f = 0;
    UInt index_filter = 0;
    UInt nn = 0;


    UInt nb_quad_facet = model->getFEEngine("FacetsFEEngine").getNbQuadraturePoints(type_facet);
    UInt nb_facet = f_filter.getSize();
    if (nb_facet == 0) continue;

    Array<Real>::const_iterator<Real> sigma_lim_it = sigma_lim.begin();

    Matrix<Real> stress_tmp(spatial_dimension, spatial_dimension);
    Matrix<Real> normal_traction(spatial_dimension, nb_quad_facet);
    Vector<Real> stress_check(nb_quad_facet);
    UInt sp2 = spatial_dimension * spatial_dimension;

    const Array<Real> & tangents = model->getTangents(type_facet);
    const Array<Real> & normals
      = model->getFEEngine("FacetsFEEngine").getNormalsOnQuadPoints(type_facet);
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
      if (stress_check.mean() > (*sigma_lim_it - tolerance)) {

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
          Real ratio = stress_check.mean()/(*sigma_lim_it);
          if (ratio > max_ratio){
            std::cout << "ratio = " << ratio << std::endl;
            ++nn;
            max_ratio = ratio;
            index_f = f;
            index_filter = f_filter(f);
          }
        }
      }
    }

    /// insertion of only 1 cohesive element in case of implicit approach. The one subjected to the highest stress.
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

    for (UInt q = 0; q < new_nb_quad_points; ++q) {
      sig_c_eff(old_nb_quad_points + q) = new_sigmas[q];
      del_c(old_nb_quad_points + q) = new_delta_c[q];
      for (UInt dim = 0; dim < spatial_dimension; ++dim) {
 	ins_stress(old_nb_quad_points + q, dim) = new_normal_traction[q](dim);
 	trac_old(old_nb_quad_points + q, dim) = new_normal_traction[q](dim);
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
  Array<Real>::iterator< Vector<Real> > traction_it =
    tractions(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::iterator< Vector<Real> > opening_it =
    opening(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::iterator< Vector<Real> > contact_traction_it =
    contact_tractions(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::iterator< Vector<Real> > contact_opening_it =
    contact_opening(el_type, ghost_type).begin(spatial_dimension);

  Array<Real>::const_iterator< Vector<Real> > normal_it =
    normal.begin(spatial_dimension);

  Array<Real>::iterator< Vector<Real> >traction_end =
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

  Array<Real>::iterator<Vector<Real> > insertion_stress_it =
    insertion_stress(el_type, ghost_type).begin(spatial_dimension);


  Real * memory_space = new Real[2*spatial_dimension];
  Vector<Real> normal_opening(memory_space, spatial_dimension);
  Vector<Real> tangential_opening(memory_space + spatial_dimension,
                                  spatial_dimension);

  /// loop on each quadrature point
  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++normal_it, ++sigma_c_it,
         ++delta_max_it, ++delta_c_it, ++damage_it, ++contact_traction_it,
         ++insertion_stress_it, ++contact_opening_it, ++delta_max_prev_it) {

    if (!model->isExplicit())
      *delta_max_it = *delta_max_prev_it;

    /// compute normal and tangential opening vectors
    Real normal_opening_norm = opening_it->dot(*normal_it);
    normal_opening  = (*normal_it);
    normal_opening *= normal_opening_norm;

    tangential_opening  = *opening_it;
    tangential_opening -=  normal_opening;

    Real tangential_opening_norm = tangential_opening.norm();

    /**
     * compute effective opening displacement
     * @f$ \delta = \sqrt{
     * \frac{\beta^2}{\kappa^2} \Delta_t^2 + \Delta_n^2 } @f$
     */
    Real delta = tangential_opening_norm * tangential_opening_norm * beta2_kappa2;

    bool penetration = normal_opening_norm < -Math::getTolerance();

    if (penetration) {
      /// use penalty coefficient in case of penetration
      *contact_traction_it = normal_opening;
      *contact_traction_it *= penalty;
      *contact_opening_it = normal_opening;
      /// don't consider penetration contribution for delta
      *opening_it = tangential_opening;
      normal_opening.clear();
    }
    else {
      delta += normal_opening_norm * normal_opening_norm;
      contact_traction_it->clear();
      contact_opening_it->clear();
    }

    delta = std::sqrt(delta);

    /// update maximum displacement and damage
    *delta_max_it = std::max(*delta_max_it, delta);
    *damage_it = std::min(*delta_max_it / *delta_c_it, 1.);

    /**
     * Compute traction @f$ \mathbf{T} = \left(
     * \frac{\beta^2}{\kappa} \Delta_t \mathbf{t} + \Delta_n
     * \mathbf{n} \right) \frac{\sigma_c}{\delta} \left( 1-
     * \frac{\delta}{\delta_c} \right)@f$
     */

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

    }
  }

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


    /// loop on each quadrature point
    for (; delta_max_it != delta_max_end;
         ++delta_max_it, ++delta_max_prev_it, ++delta_c_it) {

      if (*delta_max_prev_it == 0)
        /// elements inserted in the last step, that has not converged
        *delta_max_it = *delta_c_it / 1000;
      else
        /// elements introduced in previous steps, for which a correct
        /// value of delta_max_prev already exists
        *delta_max_it = *delta_max_prev_it;

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

  Vector<Real> normal_opening(spatial_dimension);
  Vector<Real> tangential_opening(spatial_dimension);

  Array<Real> tang_output(spatial_dimension,spatial_dimension);

  for (; tangent_it != tangent_end;
       ++tangent_it, ++normal_it, ++opening_it, ++ delta_max_it,
         ++sigma_c_it, ++delta_c_it, ++damage_it, ++contact_opening_it) {

    /// compute normal and tangential opening vectors
    *opening_it += *contact_opening_it;
    Real normal_opening_norm = opening_it->dot(*normal_it);
    normal_opening = (*normal_it);
    normal_opening *= normal_opening_norm;

    tangential_opening = *opening_it;
    tangential_opening -= normal_opening;

    Real tangential_opening_norm = tangential_opening.norm();
    bool penetration = normal_opening_norm < -Math::getTolerance();
    Real derivative = 0;
    Real t = 0;

    Real delta = tangential_opening_norm * tangential_opening_norm * beta2_kappa2;
    delta += normal_opening_norm * normal_opening_norm;
    delta = std::sqrt(delta);

    /// Delta has to be different from 0 to have finite values of tangential stiffness.
    /// At the element insertion, delta = 0. Therefore, a fictictious value is defined,
    /// for the evaluation of the first value of K.
    if (delta < Math::getTolerance())
      delta = (*delta_c_it)/1000.;

    if (normal_opening_norm >= 0.0){
      if (delta >= *delta_max_it){
        derivative = -*sigma_c_it/(delta * delta);
        t = *sigma_c_it * (1 - delta / *delta_c_it);
      }	else if (delta < *delta_max_it){
        Real tmax = *sigma_c_it * (1 - *delta_max_it / *delta_c_it);
        t = tmax / *delta_max_it * delta;
      }
    }

    Matrix<Real> n_outer_n(spatial_dimension, spatial_dimension);
    n_outer_n.outerProduct(*normal_it, *normal_it);

    if (penetration){
      /// don't consider penetration contribution for delta
      *opening_it = tangential_opening;
      /// stiffness in compression given by the penalty parameter
      *tangent_it += n_outer_n;
      *tangent_it *= penalty;
    }

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
    *tangent_it += prov;

    /// check if the tangential stiffness matrix is symmetric
    for (UInt h = 0; h < spatial_dimension; ++h){
      for (UInt l = h; l < spatial_dimension; ++l){
        if (l > h){
          Real k_ls = (*tangent_it)[spatial_dimension*h+l];
          Real k_us =  (*tangent_it)[spatial_dimension*l+h];
          //          std::cout << "k_ls = " << k_ls << std::endl;
          //          std::cout << "k_us = " << k_us << std::endl;
          if (std::abs(k_ls) > 1e-13 && std::abs(k_us) > 1e-13){
            Real error = std::abs((k_ls - k_us) / k_us);
            if (error > 1e-13){
              std::cout << "non symmetric cohesive matrix" << std::endl;
              std::cout << "error " << error << std::endl;
            }
          }
        }
      }
    }
  }
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */


INSTANSIATE_MATERIAL(MaterialCohesiveLinear);

__END_AKANTU__
