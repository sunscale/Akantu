/**
 * @file   dumper_text.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Nov 07 2017
 *
 * @brief  implementation of text dumper
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "dumper_text.hh"
#include "communicator.hh"
#include "dumper_nodal_field.hh"
#include "mesh.hh"
#include <io_helper.hh>

namespace akantu {

/* -------------------------------------------------------------------------- */
DumperText::DumperText(const std::string & basename,
                       iohelper::TextDumpMode mode, bool parallel) {
  AKANTU_DEBUG_IN();

  this->dumper = std::make_unique<iohelper::DumperText>(mode);
  this->setBaseName(basename);

  this->setParallelContext(parallel);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DumperText::registerMesh(const Mesh & mesh,
                              __attribute__((unused)) UInt spatial_dimension,
                              __attribute__((unused)) GhostType ghost_type,
                              __attribute__((unused))
                              ElementKind element_kind) {

  registerField("position",
                std::make_shared<dumpers::NodalField<Real>>(mesh.getNodes()));

  // in parallel we need node type
  UInt nb_proc = mesh.getCommunicator().getNbProc();
  if (nb_proc > 1) {
    registerField("nodes_type", std::make_shared<dumpers::NodalField<NodeFlag>>(
                                    mesh.getNodesFlags()));
  }
}

/* -------------------------------------------------------------------------- */
void DumperText::registerFilteredMesh(
    const Mesh & mesh,
    __attribute__((unused)) const ElementTypeMapArray<UInt> & elements_filter,
    const Array<UInt> & nodes_filter,
    __attribute__((unused)) UInt spatial_dimension,
    __attribute__((unused)) GhostType ghost_type,
    __attribute__((unused)) ElementKind element_kind) {

  registerField("position", std::make_shared<dumpers::NodalField<Real, true>>(
                                mesh.getNodes(), 0, 0, &nodes_filter));

  // in parallel we need node type
  UInt nb_proc = mesh.getCommunicator().getNbProc();
  if (nb_proc > 1) {
    registerField("nodes_type",
                  std::make_shared<dumpers::NodalField<NodeFlag, true>>(
                      mesh.getNodesFlags(), 0, 0, &nodes_filter));
  }
}

/* -------------------------------------------------------------------------- */
void DumperText::setBaseName(const std::string & basename) {
  AKANTU_DEBUG_IN();

  DumperIOHelper::setBaseName(basename);
  static_cast<iohelper::DumperText *>(this->dumper.get())
      ->setDataSubDirectory(this->filename + "-DataFiles");
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DumperText::setPrecision(UInt prec) {
  AKANTU_DEBUG_IN();

  static_cast<iohelper::DumperText *>(this->dumper.get())->setPrecision(prec);

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
