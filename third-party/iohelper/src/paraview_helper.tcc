#include <cassert>

#if defined(__INTEL_COMPILER)
#pragma warning ( push )
/// remark #981: operands are evaluated in unspecified order
#pragma warning ( disable : 981 )
#endif //defined(__INTEL_COMPILER)


inline std::string ParaviewHelper::dataTypeToStr(DataType data_type) {
  std::string str;
  switch(data_type) {
  case _bool   : str = "UInt8"  ; break;
  case _uint   : str = "UInt32" ; break;
  case _int    : str = "Int32"  ; break;
  case _float  : str = "Float32"; break;
  case _double : str = "Float64"; break;
  case _uint64 : str = "UInt64" ; break;
  case _int64  : str = "Int64"  ; break;
  case _uint8  : str = "UInt8"  ; break;
  }
  return str;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void ParaviewHelper::visitField(T & visited){
  this->position_flag = false;
  switch (current_stage){
  case _s_writeFieldProperty:   writeFieldProperty(visited); break;
  case _s_writePosition:        this->position_flag = true;  /* FALLTHRU */
  // [[fallthrough]] un-comment when compiler gets it
  case _s_writeField:           writeField(visited);         break;
  case _s_writeConnectivity:    writeConnectivity(visited);  break;
  case _s_writeElemType:        writeElemType(visited);      break;
  case _s_buildOffsets:         writeOffsets(visited);       break;
  default:
    std::stringstream sstr;
    sstr << "the stage " << current_stage << " is not a known paraviewhelper stage";
    IOHELPER_THROW(sstr.str(),
		   IOHelperException::_et_unknown_visitor_stage);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void ParaviewHelper::writeFieldProperty(T & data){
  if (not static_cast<bool>(data.isHomogeneous()))
    IOHELPER_THROW(std::string("try to write field property of a non homogeneous field"),
		   IOHelperException::_et_non_homogeneous_data);

  UInt dim = data.getDim();
  std::string name = data.getName();
  PDataArray(name, dim, dataTypeToStr(data.getDataType()));
}


/* -------------------------------------------------------------------------- */
template <typename T>
inline void ParaviewHelper::writeField(T & data){
  typename T::iterator it = data.begin();
  typename T::iterator end = data.end();

  compteur = 0;
  UInt dim;
  if(data.isHomogeneous()) {
    dim = data.getDim();
    if (position_flag) {
      dim = 3;
    }
    for (; it != end; ++it) {
      pushData((*it), dim);
    }
  }
  else {
    for (; it != end; ++it) {
      pushData((*it));
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void ParaviewHelper::writeConnectivity(T & data) {
  typename T::iterator it = data.begin();
  typename T::iterator end = data.end();

  for (; it != end; ++it) {
    auto type = (ElemType)it.element_type();
    //typename T::iterator::type & n = *it;
    UInt dim = (*it).size();
    assert(nb_node_per_elem[type] == dim);
    UInt * reorder = this->write_reorder[type];
    for (UInt i = 0; i < dim; ++i) {
      this->pushDatum((*it)[reorder[i]], dim);
    }
  }
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <typename T>
inline void ParaviewHelper::writeElemType(T & data){
  typename T::iterator it = data.begin();
  typename T::iterator end = data.end();

  for (; it != end; ++it) {
    auto type = (ElemType)it.element_type();
    this->pushDatum(this->paraview_code_type[type], 1);
  }
}


/* -------------------------------------------------------------------------- */
template <typename T>
inline void ParaviewHelper::writeOffsets(T & data){

  typename T::iterator it = data.begin();
  typename T::iterator end = data.end();

  UInt count = 0;
  for (; it != end; ++it) {
    count += (*it).size();
    pushDatum(count);
  }
}

/* -------------------------------------------------------------------------- */
template <template<typename T> class Cont, typename T>
inline void ParaviewHelper::pushData(const Cont<T> & n, UInt dim){
  for (UInt i = 0; i < n.size(); ++i) {
    pushDatum<T>(n[i], dim);
  }

  for (UInt i = n.size(); i < dim; ++i) { T t = T(); this->pushDatum<T>(t, dim); }
}

/* -------------------------------------------------------------------------- */
template <template<typename T> class Cont, typename T>
inline void ParaviewHelper::pushData(const Cont<T> & n) {
  for (UInt i = 0; i < n.size(); ++i) {
    pushDatum<T>(n[i], n.size());
  }
}


/* -------------------------------------------------------------------------- */
inline void ParaviewHelper::setMode(int mode){
  bflag = BASE64 & mode;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void ParaviewHelper::writePVTU(T & per_node_data, T & per_elem_data,
				      const std::vector<std::string> & vtus){

  current_stage = _s_writeFieldProperty;

  file << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" " << std::endl;
  file << "byte_order=\"LittleEndian\">" << std::endl;
  file << " <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;

  file << "  <PPoints>" << std::endl;
  file << "   <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"";
  if (bflag == BASE64) {
    file << "binary";
  } else {
    file << "ascii";
  }
  file << "\" />" << std::endl;
  file << "  </PPoints>" << std::endl;

  file << "  <PPointData>" << std::endl;
  auto itNodeField = per_node_data.begin();
  auto endNodeField = per_node_data.end();
  for (; itNodeField != endNodeField; ++itNodeField) {
    if ((*itNodeField).first != "positions") {
      (*itNodeField).second->accept(*this);
    }
  }

  file << "  </PPointData>" << std::endl;

  file << "  <PCellData>" << std::endl;
  auto itElemField = per_elem_data.begin();
  auto endElemField = per_elem_data.end();
  for (; itElemField != endElemField ; ++itElemField) {
    std::string name = (*itElemField).first;
    if (name == "connectivities" || name == "element_type") {
      continue;
    }

    (*itElemField).second->accept(*this);
  }
  file << "  </PCellData>" << std::endl;

  for (UInt l = 0; l < vtus.size(); ++l) {
    file << "  <Piece Source=\"" << vtus[l] << "\" />" << std::endl;
  }

  file << " </PUnstructuredGrid>" << std::endl;
  file << "</VTKFile>" << std::endl;
  file.close();
}

/* -------------------------------------------------------------------------- */
template <typename T>
void ParaviewHelper::pushDataFields(T & per_node_data, T & per_elem_data){

  startPointDataList();
  {
    auto itNodeField = per_node_data.begin();
    auto endNodeField = per_node_data.end();
    for (; itNodeField != endNodeField; ++itNodeField) {
      std::string name = (*itNodeField).first;
      if (name == "positions") {
        continue;
      }
      FieldInterface & f = *(*itNodeField).second;
      startData(f.getName(), f.getDim(), dataTypeToStr(f.getDataType()));
      pushField(f);
      endData();
  }
  }
  endPointDataList();

  startCellDataList();
  {
    auto itElemField = per_elem_data.begin();
    auto endElemField = per_elem_data.end();
    for (; itElemField != endElemField; ++itElemField) {
      std::string name = (*itElemField).first;
      if (name == "connectivities" || name == "element_type") {
        continue;
      }

      FieldInterface & f = *(*itElemField).second;
      startData(f.getName(), f.getDim(), dataTypeToStr(f.getDataType()));
      pushField(f);
      endData();
  }
  }
  endCellDataList();
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline  void ParaviewHelper::pushDatum(const T & n,
				       __attribute__((unused)) UInt size){
  if (bflag == BASE64) {
    b64.push<T>(n);
  } else {
    if (compteur == 0) {
      file << "      ";
    }
    ++compteur;
    file << n << " ";
  }
}

/* -------------------------------------------------------------------------- */
template <>
inline  void ParaviewHelper::pushDatum<double>(const double & n,
					       UInt size){
  if (bflag == BASE64) {
    b64.push<double>(n);
  } else {
    if (compteur % size == 0) {
      file << "     ";
    }
    file << std::setw(22);
    file << std::setprecision(15);
    file << std::scientific;
    file << n;
    file << " ";
    ++compteur;
    if (compteur % size == 0) {
      file << std::endl;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <>
inline  void ParaviewHelper::pushDatum<ElemType>(const ElemType & n,
						 __attribute__((unused)) UInt size){
  UInt n_ = this->paraview_code_type[n];
  pushDatum<UInt>(n_);
}

/* -------------------------------------------------------------------------- */
inline void ParaviewHelper::pushPosition(FieldInterface & f){
  current_stage = _s_writePosition;
  f.accept(*this);
}

/* -------------------------------------------------------------------------- */
inline void ParaviewHelper::pushField(FieldInterface & f){
  current_stage = _s_writeField;
  f.accept(*this);
}

/* -------------------------------------------------------------------------- */
inline void ParaviewHelper::pushConnectivity(FieldInterface & f) {
  current_stage = _s_writeConnectivity;
  f.accept(*this);
}

/* -------------------------------------------------------------------------- */
inline void ParaviewHelper::pushElemType(FieldInterface & f) {
  current_stage = _s_writeElemType;
  f.accept(*this);
}


/* -------------------------------------------------------------------------- */
inline void ParaviewHelper::buildOffsets(FieldInterface & f){
  current_stage = _s_buildOffsets;
  f.accept(*this);
}
/* -------------------------------------------------------------------------- */

#if defined(__INTEL_COMPILER)
#pragma warning ( pop )
#endif //defined(__INTEL_COMPILER)
