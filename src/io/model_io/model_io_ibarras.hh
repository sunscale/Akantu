/**
 * @file   model_io_ibarras.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  ModelIO implementation for IBarras input files
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
#ifndef __AKANTU_MODEL_IO_IBARRAS_HH__
#define __AKANTU_MODEL_IO_IBARRAS_HH__

/* -------------------------------------------------------------------------- */
#include "model_io.hh"
#include "structural_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

class ModelIOIBarras : public ModelIO {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ModelIOIBarras(){};
  virtual ~ModelIOIBarras(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// read a model from the file
  virtual void read(const std::string & filename, Model & model);

  /// assign sets of member to an already constructed model
  virtual void assign_sets(const std::string & filename,
                           StructuralMechanicsModel & model);
};

} // akantu

#endif /* __AKANTU_MODEL_IO_IBARRAS_HH__ */
