/**
 * @file   element_class_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Implementation parent material for IGFEM
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#include "material_igfem.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialIGFEM::MaterialIGFEM(SolidMechanicsModel & model, const ID & id) :
  Material(model, id) {};

/* -------------------------------------------------------------------------- */

MaterialIGFEM::~MaterialIGFEM() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
