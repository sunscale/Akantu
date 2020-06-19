/**
 * @file   dumper_iohelper_paraview.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Sep 26 2010
 * @date last modification: Mon Jan 22 2018
 *
 * @brief  implementations of DumperParaview
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
#include "dumper_iohelper_paraview.hh"
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <io_helper.hh>
/* -------------------------------------------------------------------------- */

namespace akantu {

DumperParaview::DumperParaview(const std::string & filename,
                               const std::string & directory, bool parallel)
    : DumperIOHelper() {
  dumper = std::make_unique<iohelper::DumperParaview>();
  setBaseName(filename);

  this->setParallelContext(parallel);

  dumper->setMode(iohelper::BASE64);
  dumper->setPrefix(directory);
  dumper->init();
}

/* -------------------------------------------------------------------------- */
DumperParaview::~DumperParaview() = default;

/* -------------------------------------------------------------------------- */
void DumperParaview::setBaseName(const std::string & basename) {
  DumperIOHelper::setBaseName(basename);
  static_cast<iohelper::DumperParaview *>(dumper.get())
      ->setVTUSubDirectory(filename + "-VTU");
}

} // namespace akantu
