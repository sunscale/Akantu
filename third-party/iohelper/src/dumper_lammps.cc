/**
 * @file   dumper_lammps.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date creation: Thu Nov 25 2010
 * @date last modification: Mon Jun 10 2013
 *
 * @brief  implementation of lammps dumper
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * IOHelper is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include <sstream>
#include "iohelper_common.hh"
#include "dumper_lammps.hh"
/* -------------------------------------------------------------------------- */

#if defined(__INTEL_COMPILER)
/// remark #981: operands are evaluated in unspecified order
#pragma warning ( disable : 981 )
#endif //defined(__INTEL_COMPILER)


namespace iohelper {

template<LammpsAtomStyle style>
DumperLammps<style>::DumperLammps(Real * bounds, const std::string & prefix)
  :Dumper(prefix), bounds(bounds) {
  this->registerDumpOptions("lammps","","", _df_proc_id);
}

/* -------------------------------------------------------------------------- */

template<LammpsAtomStyle style>
void DumperLammps<style>::dump(const std::string & current_name, const UInt count) {
  Dumper::dump(current_name, count);

  std::string filename = this->getAbsoluteFilePath(this->getBaseName(), "lammps");
  std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::trunc; 
  this->lammps_dump_file.open(filename.c_str(), mode);

  this->dumpHead(this->bounds);
  this->dumpAdd();
  this->dumpFinalize();

  this->incDumpStep();
}

/* -------------------------------------------------------------------------- */
template<LammpsAtomStyle style>
void DumperLammps<style>::dumpHead(Real * bounds) {

  this->curr_nb_atom = 0;

  if (!this->lammps_dump_file.good()) {
    std::cerr << "hach" << std::endl;
    if (this->lammps_dump_file.rdstate() & std::fstream::eofbit) {
      std::cerr << " 1 " << std::endl;
    }
    if (this->lammps_dump_file.rdstate() & std::fstream::failbit) {
      std::cerr << " 2 " << std::endl;
    }
    if (this->lammps_dump_file.rdstate() & std::fstream::badbit) {
      std::cerr << " 3 " << std::endl;
    }
    if (this->lammps_dump_file.rdstate() & std::fstream::goodbit) {
      std::cerr << " 4 " << std::endl;
    }
    exit(-1);
  }

  this->lammps_dump_file << "LAMMPS data file" << std::endl << std::endl << std::endl;
  this->nb_atom_position = lammps_dump_file.tellp();
  //dump whitespaces to later fill in nb_atoms
  this->lammps_dump_file << "                                     " << std::endl;
  this->lammps_dump_file << "0 bonds" << std::endl
                         << "1 atom types" << std::endl
                         << "0 bond types" << std::endl;
  if (bounds != NULL) {
    this->lammps_dump_file << std::endl;
    this->lammps_dump_file << bounds[0] << " " << bounds[1] <<" xlo xhi" << std::endl;
    this->lammps_dump_file << bounds[2] << " " << bounds[3] <<" ylo yhi" << std::endl;
    this->lammps_dump_file << bounds[4] << " " << bounds[5] <<" zlo zhi" << std::endl;
    this->lammps_dump_file << std::endl;
  }
  this->lammps_dump_file << "Atoms" << std::endl << std::endl;
  this->lammps_dump_file.setf(std::ios::scientific, std::ios::floatfield);
  this->lammps_dump_file.precision(16);
}

/* -------------------------------------------------------------------------- */
template<LammpsAtomStyle style>
void DumperLammps<style>::dumpAdd(int _grain_id) {
  this->grain_id = _grain_id;
  
  per_node_data[std::string("positions")]->accept(*this);

  // ContainerArray<Real> cont(_points,_dimension,_nb);
  // Field<ContainerArray<Real> > field(cont,"temporary");
  // // set temporary values
  // field.accept(*this);//this brings me to visit
}


/* -------------------------------------------------------------------------- */
template<LammpsAtomStyle style>
void DumperLammps<style>::dumpFinalize(){
  this->lammps_dump_file.seekp(this->nb_atom_position);
  this->lammps_dump_file << curr_nb_atom << " atoms";
  this->lammps_dump_file.close();
}

/* -------------------------------------------------------------------------- */

template class DumperLammps<atomic>;
template class DumperLammps<bond>;

}
