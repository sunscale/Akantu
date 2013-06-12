/**
 * @file   manual_restart.cc
 * @author Dana Christen <dana.christen@epfl.ch>
 * @date   May 15, 2013
 */

/* -------------------------------------------------------------------------- */
#include "manual_restart.hh"

//#include <iostream>
#include <fstream>

using namespace akantu;

void dumpArray(const Array<Real> & array, const std::string & fname) {
  std::ofstream outFile;
  outFile.open(fname.c_str());
  outFile.precision(9);
  outFile.setf(std::ios::scientific);
  UInt size = array.getSize();
  UInt nb_component = array.getNbComponent();
  outFile << size << std::endl;
  outFile << nb_component << std::endl;
  Array<Real>::const_iterator< Vector<Real> > tit  = array.begin(nb_component);
  Array<Real>::const_iterator< Vector<Real> > tend = array.end(nb_component);
  for(; tit != tend; ++tit) {
    for (UInt c=0; c<nb_component; ++c) {
      if (c != 0) 
	outFile << " ";
      outFile << (*tit)(c);
    }
    outFile << std::endl;
  }
  outFile.close();
}

void loadArray(Array<Real> & array, const std::string & fname) {
  std::ifstream inFile;
  inFile.open(fname.c_str());
  inFile.precision(9);
  inFile.setf(std::ios::scientific);
  UInt size(0), nb_comp(0);
  inFile >> size;
  inFile >> nb_comp;
  AKANTU_DEBUG_ASSERT(array.getNbComponent() == nb_comp, "BAD NUM OF COMPONENTS");
  Array<Real>::iterator< Vector<Real> > tit = array.begin(nb_comp);
  Array<Real>::iterator< Vector<Real> > tend = array.end(nb_comp);
  array.resize(size);
  for(UInt i(0); i < size; ++i, ++tit) {
    for (UInt c=0; c < nb_comp; ++c) {
      inFile >> (*tit)(c);
    }
  }
  inFile.close();
}
