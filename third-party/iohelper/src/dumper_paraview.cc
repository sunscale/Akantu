/**
 * @file   dumper_paraview.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Mon Jun 10 2013
 *
 * @brief  implementation of the paraview dumper
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
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <iomanip>
#include "dumper_paraview.hh"
#include "file_manager.hh"
/* -------------------------------------------------------------------------- */
namespace iohelper {

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#if defined(__INTEL_COMPILER)
/// remark #981: operands are evaluated in unspecified order
#pragma warning ( disable : 981 )
/// remark #383: value copied to temporary, reference to temporary used
#pragma warning ( disable : 383 )
#endif //defined(__INTEL_COMPILER)

/* -------------------------------------------------------------------------- */
DumperParaview::DumperParaview(const std::string & prefix) :
  Dumper(prefix), flag_compressed(false)  {
  mode_offset = 0;
  this->mode = BASE64;

  this->registerDumpOptions("pvtu", "", ".pvtu", _df_counter);
  this->registerDumpOptions("vtu", "", ".vtu", _df_counter | _df_proc_id);
  this->registerDumpOptions("pvd", "", ".pvd", _df_no_flag);
}

/* -------------------------------------------------------------------------- */
DumperParaview::~DumperParaview() {
}

/* -------------------------------------------------------------------------- */
void DumperParaview::setConnectivity(int * connectivity,
				     ElemType element_type,
				     UInt nb_elem,
				     int mode) {
  Dumper::setConnectivity(connectivity,element_type,nb_elem,mode);
  if (connectivity_mode == FORTRAN_MODE) {
    mode_offset = -1;
  } else {
    mode_offset = 0;
  }
}


/* -------------------------------------------------------------------------- */
void DumperParaview::setMode(int mode){
  Dumper::setMode(mode);
}

/* -------------------------------------------------------------------------- */
void DumperParaview::dump(const std::string & current_name, UInt count) {
  Dumper::dump(current_name, count);
  File file;

  //create directory to the files
  iohelper_mkdir(std::string(this->getAbsoluteFolderPath("pvtu")).c_str(), 0755);
  //create directory to store the per proc vtus
  iohelper_mkdir(std::string(this->getAbsoluteFolderPath("vtu")).c_str(), 0755);

  if (my_rank == root_rank){
    std::string fname = this->getAbsoluteFilePath(this->getBaseName(),
						  "pvtu");
    file.open(fname);
    if (!file.is_open()) IOHELPER_THROW("could not open file " << fname,
					_et_file_error);

    ParaviewHelper paraHelper(file, this->mode);

    std::vector<std::string> vtus;
    for (Int i = 0; i < world_size; ++i) {
      vtus.push_back(this->getRelativeFilePath(this->getBaseName(), "vtu", i));
    }

    paraHelper.writePVTU(per_node_data,per_element_data, vtus);
  }

  //open current file
  std::string fullvtuname = this->getAbsoluteFilePath(this->getBaseName(), "vtu");
  file.open(fullvtuname, std::fstream::out, flag_compressed);
  if (!file.is_open()) IOHELPER_THROW("could not open file " << fullvtuname,
				      _et_file_error);
  ParaviewHelper paraHelper(file, this->mode);

  auto pos_it = per_node_data.find("positions");
  if (pos_it == per_node_data.end())
    IOHELPER_THROW("positions field was not specified",
		   _et_missing_field);
  FieldInterface & positions = *(pos_it->second);
  UInt nb_nodes = positions.size();

  //connectivity dump
  bool point_set_flag = false;

  auto conn_it = per_element_data.find("connectivities");
  if (conn_it == per_element_data.end()) {
    point_set_flag = true;
  }
    // IOHELPER_THROW("connectivities were not specified",
    // 		   _et_missing_field);
  
  FieldInterface * connectivities_ptr = NULL;
  if (!point_set_flag){
    connectivities_ptr = (conn_it->second);
    UInt nb_elems = connectivities_ptr->size();
    
    paraHelper.writeHeader(nb_nodes, nb_elems);
  } else {
    paraHelper.writeHeader(nb_nodes, 1);
  }

  FieldInterface & connectivities = *connectivities_ptr;
  

  
  paraHelper.startDofList(3);
  paraHelper.pushPosition(positions);
  paraHelper.endDofList();

  // start to push info on cells
  paraHelper.startCells();

  // push connectivities
  paraHelper.startCellsConnectivityList();
  if (point_set_flag) {
    for (UInt i = 0; i < positions.size(); ++i) {
      paraHelper.pushDatum(i);
    }
  } else {
    paraHelper.pushConnectivity(connectivities);
  }

  paraHelper.endCellsConnectivityList();

  // build offsets
  paraHelper.startCellsoffsetsList();
  if (point_set_flag) {
    paraHelper.pushDatum(positions.size());
  } else {
    paraHelper.buildOffsets(connectivities);
  }

  paraHelper.endCellsoffsetsList();

  // push cell types
  paraHelper.startCellstypesList();
  if (point_set_flag) {
    paraHelper.pushDatum(POINT_SET);
  } else {
    auto elty_it = per_element_data.find("element_type");
    if (elty_it != per_element_data.end()) {
      paraHelper.pushField(*(elty_it->second));
    } else {
      paraHelper.pushElemType(connectivities);
    }
  }

  paraHelper.endCellstypesList();

  //close cell section
  paraHelper.endCells();

  //! add data fields
  paraHelper.pushDataFields(per_node_data, per_element_data);

  paraHelper.write_conclusion();
  file.close();

  if(my_rank == root_rank) {
    if(write_time_desc_file) {
      this->pvtu_file_names.push_back(std::make_pair(current_time,
                                                     this->getRelativeFilePath(this->getBaseName(),
                                                                               "pvtu")));
      ParaviewHelper::writeTimePVD(this->getAbsoluteFilePath(this->time_description_file_name, "pvd"),
				   this->pvtu_file_names);

      current_time += time_step;
    }
  }

  this->incDumpStep();
}

/* -------------------------------------------------------------------------- */
void DumperParaview::setVTUSubDirectory(const std::string & sub){
  this->getDumpOptions("vtu").setFolder(sub);
}

/* -------------------------------------------------------------------------- */

}
