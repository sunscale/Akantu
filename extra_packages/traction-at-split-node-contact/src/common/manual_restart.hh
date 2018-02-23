/**
 * @file   manual_restart.hh
 *
 *
 *
 * @brief
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/**
 * @file   manual_restart.hh
 * @author Dana Christen <dana.christen@epfl.ch>
 * @date   May 15, 2013
 */

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "solid_mechanics_model.hh"

void dumpArray(const akantu::Array<akantu::Real> & array,
               const std::string & fname);

void loadArray(akantu::Array<akantu::Real> & array, const std::string & fname);
void loadRestart(akantu::SolidMechanicsModel & model, const std::string & fname,
                 akantu::UInt prank);
void loadRestart(akantu::SolidMechanicsModel & model,
                 const std::string & fname);
void dumpRestart(akantu::SolidMechanicsModel & model, const std::string & fname,
                 akantu::UInt prank);
void dumpRestart(akantu::SolidMechanicsModel & model,
                 const std::string & fname);
