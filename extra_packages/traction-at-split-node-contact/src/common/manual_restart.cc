/**
 * @file   manual_restart.cc
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
  Array<Real>::const_iterator<Vector<Real>> tit = array.begin(nb_component);
  Array<Real>::const_iterator<Vector<Real>> tend = array.end(nb_component);
  for (; tit != tend; ++tit) {
    for (UInt c = 0; c < nb_component; ++c) {
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
  AKANTU_DEBUG_ASSERT(array.getNbComponent() == nb_comp,
                      "BAD NUM OF COMPONENTS");
  AKANTU_DEBUG_ASSERT(array.getSize() == size,
                      "loadArray: number of data points in file ("
                          << size << ") does not correspond to array size ("
                          << array.getSize() << ")!!");
  Array<Real>::iterator<Vector<Real>> tit = array.begin(nb_comp);
  Array<Real>::iterator<Vector<Real>> tend = array.end(nb_comp);
  array.resize(size);
  for (UInt i(0); i < size; ++i, ++tit) {
    for (UInt c = 0; c < nb_comp; ++c) {
      inFile >> (*tit)(c);
    }
  }
  inFile.close();
}

/* -------------------------------------------------------------------------- */
void loadRestart(akantu::SolidMechanicsModel & model, const std::string & fname,
                 akantu::UInt prank) {

  const akantu::Mesh & mesh = model.getMesh();
  const akantu::UInt spatial_dimension = model.getMesh().getSpatialDimension();

  const_cast<DOFSynchronizer &>(model.getDOFSynchronizer())
      .initScatterGatherCommunicationScheme();

  if (prank == 0) {
    akantu::Array<akantu::Real> full_reload_array(mesh.getNbGlobalNodes(),
                                                  spatial_dimension);
    loadArray(full_reload_array, fname);
    model.getDOFSynchronizer().scatter(model.getDisplacement(), 0,
                                       &full_reload_array);
  } else {
    model.getDOFSynchronizer().scatter(model.getDisplacement(), 0);
  }
}

/* -------------------------------------------------------------------------- */
void loadRestart(akantu::SolidMechanicsModel & model,
                 const std::string & fname) {
  loadArray(model.getDisplacement(), fname);
}

/* -------------------------------------------------------------------------- */
void dumpRestart(akantu::SolidMechanicsModel & model, const std::string & fname,
                 akantu::UInt prank) {

  const akantu::Mesh & mesh = model.getMesh();
  const akantu::UInt spatial_dimension = model.getMesh().getSpatialDimension();

  const_cast<DOFSynchronizer &>(model.getDOFSynchronizer())
      .initScatterGatherCommunicationScheme();

  if (prank == 0) {
    akantu::Array<akantu::Real> full_array(mesh.getNbGlobalNodes(),
                                           spatial_dimension);
    model.getDOFSynchronizer().gather(model.getDisplacement(), 0, &full_array);
    dumpArray(full_array, fname);
  } else {
    model.getDOFSynchronizer().gather(model.getDisplacement(), 0);
  }
}

/* -------------------------------------------------------------------------- */
void dumpRestart(akantu::SolidMechanicsModel & model,
                 const std::string & fname) {
  dumpArray(model.getDisplacement(), fname);
}
