/**
 * @file   solid_mechanics_model_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  Implementation of the inline functions of the SolidMechanicsModel class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
inline Material & SolidMechanicsModel::getMaterial(UInt mat_index) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(mat_index < materials.size(),
                     "The model " << id << " has no material no "<< mat_index);
  AKANTU_DEBUG_OUT();
  return *materials[mat_index];
}

/* -------------------------------------------------------------------------- */
inline const Material & SolidMechanicsModel::getMaterial(UInt mat_index) const {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(mat_index < materials.size(),
                     "The model " << id << " has no material no "<< mat_index);
  AKANTU_DEBUG_OUT();
  return *materials[mat_index];
}

/* -------------------------------------------------------------------------- */
inline Material & SolidMechanicsModel::getMaterial(const std::string & name) {
  AKANTU_DEBUG_IN();
  std::map<std::string, UInt>::const_iterator it = materials_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != materials_names_to_id.end(),
                     "The model " << id << " has no material named "<< name);
  AKANTU_DEBUG_OUT();
  return *materials[it->second];
}

/* -------------------------------------------------------------------------- */
inline UInt SolidMechanicsModel::getMaterialIndex(const std::string & name) const {
  AKANTU_DEBUG_IN();
  std::map<std::string, UInt>::const_iterator it = materials_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != materials_names_to_id.end(),
                      "The model " << id << " has no material named "<< name);
  AKANTU_DEBUG_OUT();
  return it->second;
}

/* -------------------------------------------------------------------------- */
inline const Material & SolidMechanicsModel::getMaterial(const std::string & name) const {
  AKANTU_DEBUG_IN();
  std::map<std::string, UInt>::const_iterator it = materials_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != materials_names_to_id.end(),
                      "The model " << id << " has no material named "<< name);
  AKANTU_DEBUG_OUT();
  return *materials[it->second];
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::setMaterialSelector(MaterialSelector & selector) {
  if(is_default_material_selector) delete material_selector;
  material_selector = &selector;
  is_default_material_selector = false;
}

/* -------------------------------------------------------------------------- */
inline FEEngine & SolidMechanicsModel::getFEEngineBoundary(const ID & name) {
  return dynamic_cast<FEEngine &>(getFEEngineClassBoundary<MyFEEngineType>(name));
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::splitElementByMaterial(const Array<Element> & elements,
                                                       Array<Element> * elements_per_mat) const {
  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  const Array<UInt> * elem_mat = NULL;

  Array<Element>::const_iterator<Element> it  = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    Element el = *it;

    if(el.type != current_element_type || el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type   = el.ghost_type;
      elem_mat = &element_index_by_material(el.type, el.ghost_type);
    }

    UInt old_id = el.element;
    el.element = (*elem_mat)(old_id, 1);
    elements_per_mat[(*elem_mat)(old_id, 0)].push_back(el);
  }
}

/* -------------------------------------------------------------------------- */
inline UInt SolidMechanicsModel::getNbDataForElements(const Array<Element> & elements,
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

  switch(tag) {
  case _gst_material_id: {
    size += elements.getSize() * 2 * sizeof(UInt);
    break;
  }
  case _gst_smm_mass: {
    size += nb_nodes_per_element * sizeof(Real) * spatial_dimension; // mass vector
    break;
  }
  case _gst_smm_for_gradu: {
    size += nb_nodes_per_element * spatial_dimension * sizeof(Real); // displacement
   break;
  }
  case _gst_smm_boundary: {
    // force, displacement, boundary
    size += nb_nodes_per_element * spatial_dimension * (2 * sizeof(Real) + sizeof(bool));
    break;
  }
  default: {  }
  }

  if(tag != _gst_material_id) {
    Array<Element> * elements_per_mat = new Array<Element>[materials.size()];
    this->splitElementByMaterial(elements, elements_per_mat);

    for (UInt i = 0; i < materials.size(); ++i) {
      size += materials[i]->getNbDataForElements(elements_per_mat[i], tag);
    }
    delete [] elements_per_mat;
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::packElementData(CommunicationBuffer & buffer,
                                                const Array<Element> & elements,
                                                SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  switch(tag) {
  case _gst_material_id: {
    packElementalDataHelper(element_index_by_material, buffer, elements, false, getFEEngine());
    break;
  }
  case _gst_smm_mass: {
    packNodalDataHelper(*mass, buffer, elements, mesh);
    break;
  }
  case _gst_smm_for_gradu: {
    packNodalDataHelper(*displacement, buffer, elements, mesh);
    break;
  }
  case _gst_smm_boundary: {
    packNodalDataHelper(*force, buffer, elements, mesh);
    packNodalDataHelper(*velocity, buffer, elements, mesh);
    packNodalDataHelper(*blocked_dofs, buffer, elements, mesh);
    break;
  }
  default: {
  }
  }

  if(tag != _gst_material_id) {
    Array<Element> * elements_per_mat = new Array<Element>[materials.size()];
    splitElementByMaterial(elements, elements_per_mat);

    for (UInt i = 0; i < materials.size(); ++i) {
      materials[i]->packElementData(buffer, elements_per_mat[i], tag);
    }

    delete [] elements_per_mat;
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::unpackElementData(CommunicationBuffer & buffer,
                                                  const Array<Element> & elements,
                                                  SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  switch(tag) {
  case _gst_material_id: {
    unpackElementalDataHelper(element_index_by_material, buffer, elements,
                             false, getFEEngine());
    break;
  }
  case _gst_smm_mass: {
    unpackNodalDataHelper(*mass, buffer, elements, mesh);
    break;
  }
  case _gst_smm_for_gradu: {
    unpackNodalDataHelper(*displacement, buffer, elements, mesh);
    break;
  }
  case _gst_smm_boundary: {
    unpackNodalDataHelper(*force, buffer, elements, mesh);
    unpackNodalDataHelper(*velocity, buffer, elements, mesh);
    unpackNodalDataHelper(*blocked_dofs, buffer, elements, mesh);
    break;
  }
  default: {
  }
  }

  if(tag != _gst_material_id) {
    Array<Element> * elements_per_mat = new Array<Element>[materials.size()];
    splitElementByMaterial(elements, elements_per_mat);

    for (UInt i = 0; i < materials.size(); ++i) {
      materials[i]->unpackElementData(buffer, elements_per_mat[i], tag);
    }

    delete [] elements_per_mat;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
inline UInt SolidMechanicsModel::getNbDataToPack(SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  //  UInt nb_nodes = mesh.getNbNodes();

  switch(tag) {
  case _gst_smm_uv: {
    size += sizeof(Real) * spatial_dimension * 2;
    break;
  }
  case _gst_smm_res: {
    size += sizeof(Real) * spatial_dimension;
    break;
  }
  case _gst_smm_mass: {
    size += sizeof(Real) * spatial_dimension;
    break;
  }
  case _gst_for_dump: {
    size += sizeof(Real) * spatial_dimension * 5;
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
inline UInt SolidMechanicsModel::getNbDataToUnpack(SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  //  UInt nb_nodes = mesh.getNbNodes();

  switch(tag) {
  case _gst_smm_uv: {
    size += sizeof(Real) * spatial_dimension * 2;
    break;
  }
  case _gst_smm_res: {
    size += sizeof(Real) * spatial_dimension;
    break;
  }
  case _gst_smm_mass: {
    size += sizeof(Real) * spatial_dimension;
    break;
  }
  case _gst_for_dump: {
    size += sizeof(Real) * spatial_dimension * 5;
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
inline void SolidMechanicsModel::packData(CommunicationBuffer & buffer,
                                         const UInt index,
                                         SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  switch(tag) {
  case _gst_smm_uv: {
    Array<Real>::const_vector_iterator it_disp = displacement->begin(spatial_dimension);
    Array<Real>::const_vector_iterator it_velo = velocity->begin(spatial_dimension);
    buffer << it_disp[index];
    buffer << it_velo[index];
    break;
  }
  case _gst_smm_res: {
    Array<Real>::const_vector_iterator it_res = residual->begin(spatial_dimension);
    buffer << it_res[index];
    break;
  }
  case _gst_smm_mass: {
    AKANTU_DEBUG_INFO("pack mass of node " << index << " which is " << (*mass)(index,0));
    Array<Real>::const_vector_iterator it_mass = mass->begin(spatial_dimension);
    buffer << it_mass[index];
    break;
  }
  case _gst_for_dump: {
    Array<Real>::const_vector_iterator it_disp = displacement->begin(spatial_dimension);
    Array<Real>::const_vector_iterator it_velo = velocity->begin(spatial_dimension);
    Array<Real>::const_vector_iterator it_acce = acceleration->begin(spatial_dimension);
    Array<Real>::const_vector_iterator it_resi = residual->begin(spatial_dimension);
    Array<Real>::const_vector_iterator it_forc = force->begin(spatial_dimension);
    buffer << it_disp[index];
    buffer << it_velo[index];
    buffer << it_acce[index];
    buffer << it_resi[index];
    buffer << it_forc[index];
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::unpackData(CommunicationBuffer & buffer,
                                           const UInt index,
                                           SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  switch(tag) {
  case _gst_smm_uv: {
    Array<Real>::vector_iterator it_disp = displacement->begin(spatial_dimension);
    Array<Real>::vector_iterator it_velo = velocity->begin(spatial_dimension);
    buffer >> it_disp[index];
    buffer >> it_velo[index];
    break;
  }
  case _gst_smm_res: {
    Array<Real>::vector_iterator it_res = residual->begin(spatial_dimension);
    buffer >> it_res[index];
    break;
  }
  case _gst_smm_mass: {
    AKANTU_DEBUG_INFO("mass of node " << index << " was " << (*mass)(index,0));
    Array<Real>::vector_iterator it_mass = mass->begin(spatial_dimension);
    buffer >> it_mass[index];
    AKANTU_DEBUG_INFO("mass of node " << index << " is now " << (*mass)(index,0));
    break;
  }
  case _gst_for_dump: {
    Array<Real>::vector_iterator it_disp = displacement->begin(spatial_dimension);
    Array<Real>::vector_iterator it_velo = velocity->begin(spatial_dimension);
    Array<Real>::vector_iterator it_acce = acceleration->begin(spatial_dimension);
    Array<Real>::vector_iterator it_resi = residual->begin(spatial_dimension);
    Array<Real>::vector_iterator it_forc = force->begin(spatial_dimension);
    buffer >> it_disp[index];
    buffer >> it_velo[index];
    buffer >> it_acce[index];
    buffer >> it_resi[index];
    buffer >> it_forc[index];
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__

#include "sparse_matrix.hh"
#include "solver.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <NewmarkBeta::IntegrationSchemeCorrectorType type>
void SolidMechanicsModel::solve(Array<Real> &increment, Real block_val,
                                bool need_factorize, bool has_profile_changed,
                                const Array<Real> &rhs) {
  
  if(has_profile_changed) {
    this->initJacobianMatrix();
    need_factorize = true;
  }

  updateResidualInternal(); //doesn't do anything for static

  if(need_factorize) {
    Real c = 0.,d = 0.,e = 0.;

    if(method == _static) {
      AKANTU_DEBUG_INFO("Solving K inc = r");
      e = 1.;
    } else {
      AKANTU_DEBUG_INFO("Solving (c M + d C + e K) inc = r");

      NewmarkBeta * nmb_int = dynamic_cast<NewmarkBeta *>(integrator);
      c = nmb_int->getAccelerationCoefficient<type>(time_step);
      d = nmb_int->getVelocityCoefficient<type>(time_step);
      e = nmb_int->getDisplacementCoefficient<type>(time_step);
    }


    jacobian_matrix->clear();
    // J = c M + d C + e K
    if(stiffness_matrix)
      jacobian_matrix->add(*stiffness_matrix, e);

    if(mass_matrix)
      jacobian_matrix->add(*mass_matrix, c);

#if !defined(AKANTU_NDEBUG)
    if(mass_matrix && AKANTU_DEBUG_TEST(dblDump))
      mass_matrix->saveMatrix("M.mtx");
#endif

    if(velocity_damping_matrix)
      jacobian_matrix->add(*velocity_damping_matrix, d);

    jacobian_matrix->applyBoundary(*blocked_dofs, block_val);

#if !defined(AKANTU_NDEBUG)
    if(AKANTU_DEBUG_TEST(dblDump))
      jacobian_matrix->saveMatrix("J.mtx");
#endif
    solver->factorize();
  }

  if (rhs.getSize() != 0)
    solver->setRHS(rhs);
  else
    solver->setRHS(*residual);

  // solve @f[ J \delta w = r @f]
  solver->solve(increment);

  UInt nb_nodes = displacement-> getSize();
  UInt nb_degree_of_freedom = displacement->getNbComponent() * nb_nodes;

  bool * blocked_dofs_val = blocked_dofs->storage();
  Real * increment_val = increment.storage();

  for (UInt j = 0; j < nb_degree_of_freedom;
       ++j,++increment_val, ++blocked_dofs_val) {
    if ((*blocked_dofs_val))
      *increment_val = 0.0;
    }

}


/* -------------------------------------------------------------------------- */
template<SolveConvergenceMethod cmethod, SolveConvergenceCriteria criteria>
bool SolidMechanicsModel::solveStatic(Real tolerance, UInt max_iteration,
				      bool do_not_factorize) {

  AKANTU_DEBUG_INFO("Solving Ku = f");
  AKANTU_DEBUG_ASSERT(stiffness_matrix != NULL,
                      "You should first initialize the implicit solver and assemble the stiffness matrix by calling initImplicit");

  AnalysisMethod analysis_method=method;
  Real error = 0.;
  method=_static;
  bool converged = this->template solveStep<cmethod, criteria>(tolerance, error, max_iteration, do_not_factorize);
  method=analysis_method;
  return converged;

}

/* -------------------------------------------------------------------------- */
template<SolveConvergenceMethod cmethod, SolveConvergenceCriteria criteria>
bool SolidMechanicsModel::solveStep(Real tolerance,
                                    UInt max_iteration) {
  Real error = 0.;
  return this->template solveStep<cmethod,criteria>(tolerance,
                                                    error,
                                                    max_iteration);
}

/* -------------------------------------------------------------------------- */
template<SolveConvergenceMethod cmethod, SolveConvergenceCriteria criteria>
bool SolidMechanicsModel::solveStep(Real tolerance, Real & error, UInt max_iteration,
				    bool do_not_factorize) {
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeSolveStepEvent(method));
  this->implicitPred();
  this->updateResidual();

  AKANTU_DEBUG_ASSERT(stiffness_matrix != NULL,
                      "You should first initialize the implicit solver and assemble the stiffness matrix");

  bool need_factorize = !do_not_factorize;

  if (method==_implicit_dynamic) {
    AKANTU_DEBUG_ASSERT(mass_matrix != NULL,
                        "You should first initialize the implicit solver and assemble the mass matrix");
  }

  switch (cmethod) {
  case _scm_newton_raphson_tangent:
  case _scm_newton_raphson_tangent_not_computed:
    break;
  case _scm_newton_raphson_tangent_modified:
    this->assembleStiffnessMatrix();
    break;
  default:
    AKANTU_DEBUG_ERROR("The resolution method " << cmethod << " has not been implemented!");
  }

  UInt iter = 0;
  bool converged = false;
  error = 0.;
  if(criteria == _scc_residual) {
    converged = this->testConvergence<criteria> (tolerance, error);
    if(converged) return converged;
  }

  do {
    if (cmethod == _scm_newton_raphson_tangent)
      this->assembleStiffnessMatrix();

    solve<NewmarkBeta::_displacement_corrector> (*increment, 1., need_factorize);

    this->implicitCorr();

    if(criteria == _scc_residual) this->updateResidual();

    converged = this->testConvergence<criteria> (tolerance, error);

    if(criteria == _scc_increment && !converged) this->updateResidual();
    //this->dump();

    iter++;
    AKANTU_DEBUG_INFO("[" << criteria << "] Convergence iteration "
                      << std::setw(std::log10(max_iteration)) << iter
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


  } while (!converged && iter < max_iteration);

  // this makes sure that you have correct strains and stresses after the solveStep function (e.g., for dumping)
  if(criteria == _scc_increment) this->updateResidual();

  if (converged) {
    EventManager::sendEvent(SolidMechanicsModelEvent::AfterSolveStepEvent(method));
  } else if(iter == max_iteration) {
    AKANTU_DEBUG_WARNING("[" << criteria << "] Convergence not reached after "
                         << std::setw(std::log10(max_iteration)) << iter <<
                         " iteration" << (iter == 1 ? "" : "s") << "!" << std::endl);
  }

  return converged;
}
