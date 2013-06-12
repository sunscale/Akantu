/**
 * @file   manual_restart.hh
 * @author Dana Christen <dana.christen@epfl.ch>
 * @date   May 15, 2013
 */

/* -------------------------------------------------------------------------- */
#include "aka_vector.hh"

void dumpArray(const akantu::Array<akantu::Real> & array, const std::string & fname);

void loadArray(akantu::Array<akantu::Real> & array, const std::string & fname);
