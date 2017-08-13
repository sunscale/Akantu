/**
 * @file   heat_transfer_model_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 *
 * @date creation: Fri Aug 20 2010
 * @date last modification: Tue Dec 08 2015
 *
 * @brief  Implementation of the inline functions of the HeatTransferModel class
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
inline FEEngine & HeatTransferModel::getFEEngineBoundary(const std::string & name) {
  return dynamic_cast<FEEngine &>(getFEEngineClassBoundary<MyFEEngineType>(name));
}

/* -------------------------------------------------------------------------- */
inline UInt HeatTransferModel::getNbDataToPack(SynchronizationTag tag) const{
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes = getFEEngine().getMesh().getNbNodes();

  switch(tag) {
  case _gst_htm_temperature:
  case _gst_htm_capacity: {
    size += nb_nodes * sizeof(Real);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline UInt HeatTransferModel::getNbDataToUnpack(SynchronizationTag tag) const{
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes = getFEEngine().getMesh().getNbNodes();

  switch(tag) {
  case _gst_htm_capacity:
  case _gst_htm_temperature: {
    size += nb_nodes * sizeof(Real);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::packData(CommunicationBuffer & buffer,
					const UInt index,
					SynchronizationTag tag) const{
  AKANTU_DEBUG_IN();

  switch(tag) {
  case _gst_htm_capacity:
    buffer << (*capacity_lumped)(index);
    break;
  case _gst_htm_temperature: {
    buffer << (*temperature)(index);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::unpackData(CommunicationBuffer & buffer,
					  const UInt index,
					  SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  switch(tag) {
  case _gst_htm_capacity: {
    buffer >> (*capacity_lumped)(index);
    break;
  }
  case _gst_htm_temperature: {
    buffer >> (*temperature)(index);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline UInt HeatTransferModel::getNbDataForElements(const Array<Element> & elements,
                                                    SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = 0;
  Array<Element>::const_iterator<Element> it  = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

#ifndef AKANTU_NDEBUG
  size += elements.size() * spatial_dimension * sizeof(Real); /// position of the barycenter of the element (only for check)
  //  size += spatial_dimension * nb_nodes_per_element * sizeof(Real); /// position of the nodes of the element
#endif

  switch(tag) {
  case _gst_htm_capacity: {
    size += nb_nodes_per_element * sizeof(Real); // capacity vector
    break;
  }
  case _gst_htm_temperature: {
    size += nb_nodes_per_element * sizeof(Real); // temperature
    break;
  }
  case _gst_htm_gradient_temperature: {
    size += getNbIntegrationPoints(elements) * spatial_dimension * sizeof(Real); // temperature gradient
    size += nb_nodes_per_element * sizeof(Real); // nodal temperatures
    //    size += spatial_dimension * nb_nodes_per_element * sizeof(Real); // shape derivatives
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::packElementData(CommunicationBuffer & buffer,
                                               const Array<Element> & elements,
                                               SynchronizationTag tag) const {
#ifndef AKANTU_NDEBUG
  Array<Element>::const_iterator<Element> bit  = elements.begin();
  Array<Element>::const_iterator<Element> bend = elements.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;
    Vector<Real> barycenter(spatial_dimension);
    mesh.getBarycenter(element.element, element.type, barycenter.storage(), element.ghost_type);
    buffer << barycenter;
  }
  // packNodalDataHelper(mesh.getNodes(), buffer, elements);
#endif

  switch (tag){
  case _gst_htm_capacity: {
    packNodalDataHelper(*capacity_lumped, buffer, elements, mesh);
    break;
  }
  case _gst_htm_temperature: {
    packNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  case _gst_htm_gradient_temperature: {
    packElementalDataHelper(temperature_gradient, buffer, elements, true, getFEEngine());
    packNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::unpackElementData(CommunicationBuffer & buffer,
                                                 const Array<Element> & elements,
                                                 SynchronizationTag tag) {
#ifndef AKANTU_NDEBUG
  Array<Element>::const_iterator<Element> bit  = elements.begin();
  Array<Element>::const_iterator<Element> bend = elements.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;

    Vector<Real> barycenter_loc(spatial_dimension);
    mesh.getBarycenter(element.element, element.type, barycenter_loc.storage(), element.ghost_type);

    Vector<Real> barycenter(spatial_dimension);
    buffer >> barycenter;
    Real tolerance = 1e-15;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      if(!(std::abs(barycenter(i) - barycenter_loc(i)) <= tolerance))
	AKANTU_DEBUG_ERROR("Unpacking an unknown value for the element: "
			   << element
			   << "(barycenter[" << i << "] = " << barycenter_loc(i)
			   << " and buffer[" << i << "] = " << barycenter(i) << ") - tag: " << tag);
    }
  }

  // Vector<Real> coords(spatial_dimension);
  // Real * nodes = getFEEngine().getMesh().getNodes().storage();
  // for (UInt n = 0; n < nb_nodes_per_element; ++n) {
  //   buffer >> coords;
  //   UInt offset_conn = conn[el_offset + n];
  //   Real * coords_local = nodes+spatial_dimension*offset_conn;
  //   for (UInt i = 0; i < spatial_dimension; ++i) {
  //     if(!(std::abs(coords(i) - coords_local[i]) <= tolerance))
  //       AKANTU_EXCEPTION("Unpacking to wrong node for the element : "
  //       		 << element
  //       		 << "(coords[" << i << "] = " << coords_local[i]
  //       		 << " and buffer[" << i << "] = " << coords(i) << ")");
  //   }
  // }
#endif

  switch (tag){
  case _gst_htm_capacity: {
    unpackNodalDataHelper(*capacity_lumped, buffer, elements, mesh);
    break;
  }
  case _gst_htm_temperature: {
    unpackNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  case _gst_htm_gradient_temperature: {
    unpackElementalDataHelper(temperature_gradient, buffer, elements, true, getFEEngine());
    unpackNodalDataHelper(*temperature, buffer, elements, mesh);

    // //    Real tolerance = 1e-15;
    // if (!Math::are_vector_equal(spatial_dimension,gtemp.storage(),it_gtemp[element.element].storage())){
    //   Real dist = Math::distance_3d(gtemp.storage(), it_gtemp[element.element].storage());
    //   debug::debugger.getOutputStream().precision(20);
    //   std::stringstream temperatures_str;
    //   temperatures_str.precision(20);
    //   temperatures_str << std::scientific << "temperatures are ";
    //   for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    //     UInt offset_conn = conn[el_offset + n];
    //     temperatures_str << (*temperature)(offset_conn) << " ";
    //   }
    //   Array<Real>::matrix_iterator it_shaped =
    //     const_cast<Array<Real> &>(getFEEngine().getShapesDerivatives(element.type, ghost_type))
    //     .begin(nb_nodes_per_element,spatial_dimension);


    //   AKANTU_EXCEPTION("packed gradient do not match for element " << element.element << std::endl
    //     	       << "buffer is " << gtemp << " local is " << it_gtemp[element.element]
    //     	       << " dist is " << dist << std::endl
    //     	       << temperatures_str.str() << std::endl
    //     	       << std::scientific << std::setprecision(20)
    //     	       << " distant temperatures " << temp_nodes
    //     	       << "temperature gradient size " << temperature_gradient(element.type, ghost_type).size()
    //     	       << " number of ghost elements " << getFEEngine().getMesh().getNbElement(element.type,_ghost)
    //     	       << std::scientific << std::setprecision(20)
    //     	       << " shaped " << shaped
    //     	       << std::scientific << std::setprecision(20)
    //     	       << " local shaped " << it_shaped[element.element]);
    // }
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */

template<SolveConvergenceMethod cmethod, SolveConvergenceCriteria criteria>
bool HeatTransferModel::solveStep(Real tolerance,
				  UInt max_iteration) {
  Real error = 0.;
  return this->template solveStep<cmethod,criteria>(tolerance,
                                                    error,
                                                    max_iteration);
}

/* -------------------------------------------------------------------------- */
template<SolveConvergenceMethod cmethod, SolveConvergenceCriteria criteria>
bool HeatTransferModel::solveStep(Real tolerance, Real & error,
				  UInt max_iteration,
				  bool do_not_factorize) {
  //EventManager::sendEvent(HeatTransferModelEvent::BeforeSolveStepEvent(method));
  this->implicitPred();
  this->updateResidual();

  AKANTU_DEBUG_ASSERT(conductivity_matrix != NULL,
                      "You should first initialize the implicit solver and assemble the conductivity matrix");

  bool need_factorize = !do_not_factorize;

  if (method == _implicit_dynamic) {
    AKANTU_DEBUG_ASSERT(capacity_matrix != NULL,
                        "You should first initialize the implicit solver and assemble the mass matrix");
  }

  switch (cmethod) {
  case _scm_newton_raphson_tangent:
  case _scm_newton_raphson_tangent_not_computed:
    break;
  case _scm_newton_raphson_tangent_modified:
    this->assembleConductivityMatrix();
    break;
  default:
    AKANTU_DEBUG_ERROR("The resolution method " << cmethod << " has not been implemented!");
  }

  this->n_iter = 0;
  bool converged = false;
  error = 0.;
  if(criteria == _scc_residual) {
    converged = this->testConvergence<criteria> (tolerance, error);
    if(converged) return converged;
  }

  do {
    if (cmethod == _scm_newton_raphson_tangent)
      this->assembleConductivityMatrix();

    solve<GeneralizedTrapezoidal::_temperature_corrector>(*increment, 1., need_factorize);

    this->implicitCorr();

    if(criteria == _scc_residual) this->updateResidual();

    converged = this->testConvergence<criteria> (tolerance, error);

    if(criteria == _scc_increment && !converged) this->updateResidual();
    //   this->dump();

    this->n_iter++;
    AKANTU_DEBUG_INFO("[" << criteria << "] Convergence iteration "
                      << std::setw(std::log10(max_iteration)) << this->n_iter
                      << ": error " << error << (converged ? " < " : " > ") << tolerance);

    switch (cmethod) {
    case _scm_newton_raphson_tangent:
      need_factorize = true;
      break;
    case _scm_newton_raphson_tangent_not_computed:
    case _scm_newton_raphson_tangent_modified:
      need_factorize = false;
      break;
    default:
      AKANTU_DEBUG_ERROR("The resolution method " << cmethod << " has not been implemented!");
    }


  } while (!converged && this->n_iter < max_iteration);

  // this makes sure that you have correct strains and stresses after the solveStep function (e.g., for dumping)
  if(criteria == _scc_increment) this->updateResidual();

  if (converged) {
    //EventManager::sendEvent(HeatTransferModelEvent::AfterSolveStepEvent(method));
  } else if(this->n_iter == max_iteration) {
    AKANTU_DEBUG_WARNING("[" << criteria << "] Convergence not reached after "
                         << std::setw(std::log10(max_iteration)) << this->n_iter <<
                         " iteration" << (this->n_iter == 1 ? "" : "s") << "!" << std::endl);
  }

  return converged;
}

/* -------------------------------------------------------------------------- */

template <GeneralizedTrapezoidal::IntegrationSchemeCorrectorType type>
void HeatTransferModel::solve(Array<Real> &increment, Real block_val,
			      bool need_factorize, bool has_profile_changed){

  updateResidualInternal(); //doesn't do anything for static

  if(need_factorize) {
    Real c = 0.,e = 0.;

    if(method == _static) {
      AKANTU_DEBUG_INFO("Solving K inc = r");
      e = 1.;
    } else {
      AKANTU_DEBUG_INFO("Solving (c M + e K) inc = r");

      GeneralizedTrapezoidal * trap_int = dynamic_cast<GeneralizedTrapezoidal*>(integrator);
      c = trap_int->getTemperatureRateCoefficient<type>(time_step);
      e = trap_int->getTemperatureCoefficient<type>(time_step);

      //      std::cout << "c " << c << " e " << e << std::endl;
    }


    jacobian_matrix->clear();
    // J = c M + e K
    if(conductivity_matrix)
      jacobian_matrix->add(*conductivity_matrix, e);

    if(capacity_matrix)
      jacobian_matrix->add(*capacity_matrix, c);

#if !defined(AKANTU_NDEBUG)
    //    if(capacity_matrix && AKANTU_DEBUG_TEST(dblDump))
      capacity_matrix->saveMatrix("M.mtx");
#endif

    jacobian_matrix->applyBoundary(*blocked_dofs, block_val);

#if !defined(AKANTU_NDEBUG)
    //if(AKANTU_DEBUG_TEST(dblDump))
      jacobian_matrix->saveMatrix("J.mtx");
#endif
    solver->factorize();
  }

  // if (rhs.size() != 0)
  //  solver->setRHS(rhs);
  // else

  solver->setOperators();

  solver->setRHS(*residual);

  // solve @f[ J \delta w = r @f]
  solver->solve(increment);

  UInt nb_nodes = temperature->size();
  UInt nb_degree_of_freedom = temperature->getNbComponent() * nb_nodes;

  bool * blocked_dofs_val = blocked_dofs->storage();
  Real * increment_val = increment.storage();

  for (UInt j = 0; j < nb_degree_of_freedom;
       ++j,++increment_val, ++blocked_dofs_val) {
    if ((*blocked_dofs_val))
      *increment_val = 0.0;
    }

}

/* -------------------------------------------------------------------------- */

