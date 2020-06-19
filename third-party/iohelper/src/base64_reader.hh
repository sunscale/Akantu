/**
 * @file   base64_reader.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Thu Nov 01 2012
 *
 * @brief  header for base64 reader
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

#ifndef BASE64_READER_H
#define BASE64_READER_H

/** Class that allow to push binary data
    in base64 format to any file.
    This class is mainly used by the paraview helper
    to create binary XML VTK files.
    The conversion is a 4/3 size conversion. */

class Base64Reader{

 public:

  Base64Reader(){
    InitBase64Stuff();
    n = 0;
    n_coded = 0;
    read_bytes_in_double=0;
  };

  //! push to file a double
  int PopDoubleInBase64(double & c);
  //! push to file an integer
  int PopIntegerInBase64(int & c);
  //! push to file a Byte
  int PopByteInBase64(unsigned char & c);
  //! decode 3 bytes from 4 Base64 bytes (4/3 ratio)
  int Decode(unsigned char c0,unsigned char c1,unsigned char c2,unsigned char c3,
	      unsigned char * r0,unsigned char * r1,unsigned char * r2);

  //! set the reading stream
  void SetStream(char * str,int l);

 private:

  //! initialisation process
  void InitBase64Stuff();
  //! decoding table
  char dtable[256];
  //! encoding table
  char etable[256];
  //! stage in conversion process(1,2 or 3)
  int n;
  //! stage for coded bytes
  int n_coded;
  //! used to code/decode
  unsigned char igroup[3],ogroup[4];


  //! starting stream position
  char * start_ptr;
  //! last reading position
  char * current_ptr;
  //! size of the stream
  int len;
  //! readed bytes
  unsigned char read_byte[3];
  //! readed bytes but still coded
  char read_coded_byte[4];


  //! temporary double
  double temp_double;
  //! readed bytes in double
  int read_bytes_in_double;
};



inline void Base64Reader::SetStream(char * str,int l){
  current_ptr = str;
  start_ptr = str;
  len = l;

  DUMP("n_coded " << n_coded,DBG_DETAIL);
  if (n_coded > 0){
    //je relie l'ancien stream et le nouveau

    if (n_coded == 1 && current_ptr < start_ptr+len) {
      read_coded_byte[1] = *current_ptr;
      ++n_coded;
      ++current_ptr;
    }
    DUMP("n_coded " << n_coded,DBG_DETAIL);
    if (n_coded == 2 && current_ptr < start_ptr+len) {
      read_coded_byte[2] = *current_ptr;
      ++n_coded;
      ++current_ptr;
    }
    DUMP("n_coded " << n_coded,DBG_DETAIL);
    if (n_coded == 3 && current_ptr < start_ptr+len) {
      read_coded_byte[3] = *current_ptr;
      ++n_coded;
      ++current_ptr;
    }
    DUMP("n_coded " << n_coded,DBG_DETAIL);
  }


};


inline int Base64Reader::PopIntegerInBase64(int & d){
  int read = 1;

  DUMP("pushing " << d << " ( n = " << n << " )",DBG_DETAIL);
  unsigned char * c = (unsigned char*)&d;
  for (unsigned int i = 0 ; i < sizeof(int) ; ++i){
    read &= PopByteInBase64(c[i]);
    if (!read) break;
  }
  return read;
}

/* inline void Base64Reader::PopStrInBase64(char * str){ */

/*   FATAL("unimplemented for the moment"); */
/* } */

inline int Base64Reader::PopDoubleInBase64(double & d){
  int read = 1;
  DUMP("pop double , already red  = " << read_bytes_in_double,DBG_DETAIL);

  unsigned char * c = (unsigned char*)&temp_double;
  for (unsigned int i = read_bytes_in_double ; i < sizeof(double) ; ++i){
    read &= PopByteInBase64(c[i]);
    if (read) ++read_bytes_in_double;
    if (!read) break;
  }

  if (read) {
    read_bytes_in_double = 0;
    d = temp_double;
    DUMP("readed double " << d,DBG_DETAIL);
  }
  return read;
}

inline int Base64Reader::PopByteInBase64(unsigned char & c){
  //initialise les blocs
  if (n < 0 || n > 3) FATAL("this should not append check source code");

  DUMP("pop byte",DBG_DETAIL);

  // si je suis a n = 3 je dois relire 4 octets de base64
  if (n == 0){

    DUMP("n_coded " << n_coded,DBG_DETAIL);
    if (n_coded == 0){
      if (current_ptr < start_ptr+len) {read_coded_byte[0] = current_ptr[0];++n_coded;DUMP("n_coded " << n_coded,DBG_DETAIL);}
      else {DUMP("lala" << n_coded,DBG_DETAIL);return 0;}
      if (current_ptr+1 < start_ptr+len) {read_coded_byte[1] = current_ptr[1];++n_coded;DUMP("n_coded " << n_coded,DBG_DETAIL);}
      else {DUMP("lala" << n_coded,DBG_DETAIL);return 0;}
      if (current_ptr+2 < start_ptr+len) {read_coded_byte[2] = current_ptr[2];++n_coded;DUMP("n_coded " << n_coded,DBG_DETAIL);}
      else {DUMP("lala" << n_coded,DBG_DETAIL);return 0;}
      if (current_ptr+3 < start_ptr+len) {read_coded_byte[3] = current_ptr[3];++n_coded;DUMP("n_coded " << n_coded,DBG_DETAIL);}
      else {DUMP("lala" << n_coded,DBG_DETAIL);return 0;}
      current_ptr+=4;
    }
  }

  if (n_coded == 4)
    Decode(read_coded_byte[0],read_coded_byte[1],read_coded_byte[2],read_coded_byte[3],
	   &read_byte[0],&read_byte[1],&read_byte[2]);

  c = read_byte[n];

  DUMP("readed charater "<< (int)c,DBG_DETAIL);
  ++n;

  if (n == 3) {
    n = 0;
    n_coded = 0;
  }
  return 1;

}


inline int Base64Reader::Decode(unsigned char c0,unsigned char c1,unsigned char c2,unsigned char c3,
				 unsigned char * r0,unsigned char * r1,unsigned char * r2){

  unsigned char d0, d1, d2, d3;

  d0 = dtable[0+c0];
  d1 = dtable[0+c1];
  d2 = dtable[0+c2];
  d3 = dtable[0+c3];



  DUMP("d0 " << (int)d0 << " d1 " << (int)d1 << " d2 " << (int)d2 << " d3 " << (int)d3,DBG_DETAIL);

  // Decode the 3 bytes

  *r0 = ((d0 << 2) & 0xFC) | ((d1 >> 4) & 0x03);
  *r1 = ((d1 << 4) & 0xF0) | ((d2 >> 2) & 0x0F);
  *r2 = ((d2 << 6) & 0xC0) | ((d3 >> 0) & 0x3F);

  DUMP("r0 " << (int)*r0 << " r1 " << (int)*r1 << " r2 " << (int)*r2,DBG_DETAIL);

  // Return the number of bytes actually decoded

  if (c2 == '=')
    {
    return 1;
    }
  if (c3 == '=')
    {
    return 2;
    }
  return 3;
}



inline void Base64Reader::InitBase64Stuff(){
  memset(dtable,0xFF,256);
  memset(etable,0xFF,256);

  for(int i=0;i<9;i++){
    etable[i]= 'A'+i;
    dtable[0+etable[i]] = i;
    etable[i+9]= 'J'+i;
    dtable[0+etable[i+9]] = i+9;
    etable[26+i]= 'a'+i;
    dtable[0+etable[26+i]] = i + 26;
    etable[26+i+9]= 'j'+i;
    dtable[0+etable[26+i+9]] = i + 26 + 9;
  }
  for(int i= 0;i<8;i++){
    etable[i+18]= 'S'+i;
    dtable[0+etable[i+18]] = i + 18;
    etable[26+i+18]= 's'+i;
    dtable[0+etable[26+i+18]] = 26 + i + 18;
  }
  for(int i= 0;i<10;i++){
    etable[52+i]= '0'+i;
    dtable[0+etable[i+52]] = i + 52;
  }
  etable[62]= '+';
  dtable[0+etable[62]] = 62;
  etable[63]= '/';
  dtable[0+etable[63]] = 63;
}

#endif //BASE64_WRITER_H

