/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateAcceleration() {
  AKANTU_DEBUG_IN();

  updateResidualInternal();

  if(method == _explicit_lumped_mass) {
    /* residual = residual_{n+1} - M * acceleration_n therefore
       solution = increment acceleration not acceleration */
    solveLumped(*increment_acceleration,
		*mass,
		*residual,
		*blocked_dofs,
		f_m2a);
  } else if (method == _explicit_consistent_mass) {
    solve<NewmarkBeta::_acceleration_corrector>(*increment_acceleration);
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::explicitPred() {
  AKANTU_DEBUG_IN();

  if(increment_flag) {
    if(previous_displacement) increment->copy(*previous_displacement);
    else increment->copy(*displacement);
  }

  AKANTU_DEBUG_ASSERT(integrator,"itegrator should have been allocated: "
		      << "have called initExplicit ? "
		      << "or initImplicit ?");

  integrator->integrationSchemePred(time_step,
				    *displacement,
				    *velocity,
				    *acceleration,
				    *blocked_dofs);

  if(increment_flag) {
    Real * inc_val = increment->storage();
    Real * dis_val = displacement->storage();
    UInt nb_degree_of_freedom = displacement->getSize() * displacement->getNbComponent();

    for (UInt n = 0; n < nb_degree_of_freedom; ++n) {
      *inc_val = *dis_val - *inc_val;
      inc_val++;
      dis_val++;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::explicitCorr() {
  AKANTU_DEBUG_IN();

  integrator->integrationSchemeCorrAccel(time_step,
					 *displacement,
					 *velocity,
					 *acceleration,
					 *blocked_dofs,
					 *increment_acceleration);

  if(previous_displacement) previous_displacement->copy(*displacement);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::solveStep() {
  AKANTU_DEBUG_IN();

  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeSolveStepEvent(method));

  this->explicitPred();
  this->updateResidual();
  this->updateAcceleration();
  this->explicitCorr();

  EventManager::sendEvent(SolidMechanicsModelEvent::AfterSolveStepEvent(method));

  AKANTU_DEBUG_OUT();
}
