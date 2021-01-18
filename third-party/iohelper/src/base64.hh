/**
 * @file   base64.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Wed Nov 13 2013
 *
 * @brief  header for base64 handling
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
#ifndef IOHELPER_BASE64_H_
#define IOHELPER_BASE64_H_
/* -------------------------------------------------------------------------- */
#include <vector>
#include "file_manager.hh"
#ifdef USING_ZLIB
#include <zlib.h>
#endif
/* -------------------------------------------------------------------------- */

namespace iohelper {

#if defined(__INTEL_COMPILER)
#pragma warning ( push )
/// remark #981: operands are evaluated in unspecified order
#pragma warning ( disable : 981 )
#endif //defined(__INTEL_COMPILER)

/** Class that allow to push binary data
    in base64 format to any file.
    This class is mainly used by the paraview helper
    to create binary XML VTK files.
    The conversion is a 4/3 size conversion. */


class Base64Writer{

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

 public:

  Base64Writer(File & f);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  //! for any packet in base64 a little header is used that is an int = nbBytes written
  void WriteHeader();
  //! this is used to allocate the memory for the final count of bytes
  void CreateHeader();
  //! when all stream is ready buffer is sent to file
  inline void DumpToFile();
  //! empty temporary buffer
  void ClearBuffer();

  template <typename T> inline void push(T t);

  //! notify that we don't want to add any data. Closing the current buffer
  void finish();
  //! decode 3 bytes from 4 Base64 bytes (4/3 ratio)
  int Decode(char c0,char c1,char c2,char c3,
	      char * r0,char * r1,char * r2);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
 private:
  //! push to file a Byte
  inline void PushByteInBase64(unsigned char c);

  //! when 4 bytes are ready they are dumped to buffer by this function
  inline void dumpToBuffer();
  //! initialisation process
  inline void InitBase64Stuff();
  //! primitive function to push bytes to file
  //void ochar(int c);

  //! decoding table
  char dtable[256];
  //! encoding table
  char etable[256];
  //! stage in conversion process(1,2 or 3)
  int n;
  //! used to code/decode
  unsigned char igroup[3],ogroup[4];

  //!unused
  int linelength;
  //!unused
  //  int maxlinelength;

  //! LM file descriptor
  File & file;

  //! buffer to cache data
  std::vector<unsigned char> buffer;
  //! number of bytes written to buffer
  long nbBytes;
  //! for rewind need o floating index
  int start;

};

/* -------------------------------------------------------------------------- */
template<typename T> inline void Base64Writer::push(T t) {
  auto * c = (unsigned char *)&t;
  for (unsigned int i = 0 ; i < sizeof(T) ; ++i){
    this->PushByteInBase64(c[i]);
  }
}

template<> inline void Base64Writer::push<char *>(char * t) {
  auto * c = (unsigned char *)t;
  for (unsigned int i = 0 ; i < 512 ; ++i){
    if (t[i] == '\0') {
      break;
    }
    PushByteInBase64(c[i]);
  }
}

inline void Base64Writer::InitBase64Stuff(){
  memset(dtable,0xFF,256);
  memset(etable,0xFF,256);

  for(int i=0;i<9;i++){
    etable[i]= (char)('A'+i);
    dtable[0+etable[i]] = (char)i;
    etable[i+9]= (char)('J'+i);
    dtable[(0+etable[i+9])] = (char)(i+9);
    etable[26+i]= (char)('a'+i);
    dtable[(0+etable[26+i])] = (char)(i + 26);
    etable[(26+i+9)]= (char)('j'+i);
    dtable[(0+etable[26+i+9])] = (char)(i + 26 + 9);
  }
  for(int i= 0;i<8;i++){
    etable[i+18]= (char)('S'+i);
    dtable[(0+etable[i+18])] = (char)(i + 18);
    etable[26+i+18]= (char)('s'+i);
    dtable[(0+etable[26+i+18])] = (char)(26 + i + 18);
  }
  for(char i= 0;i<10;i++){
    etable[52+i]= (char)('0'+i);
    dtable[(0+etable[i+52])] = (char)(i + 52);
  }
  etable[62]= '+';
  dtable[0+etable[62]] = 62;
  etable[63]= '/';
  dtable[0+etable[63]] = 63;
}
/* -------------------------------------------------------------------------- */


// inline void Base64Writer::PushIntegerInBase64(int d){
//   //  DUMP("pushing " << d << " ( n = " << n << " )",DBG_ALL);
//   unsigned char * c = (unsigned char*)&d;
//   for (unsigned int i = 0 ; i < sizeof(int) ; ++i){
//     PushByteInBase64(c[i]);
//   }
// }
// /* -------------------------------------------------------------------------- */


// inline void Base64Writer::PushStrInBase64(char * str){
//   unsigned char * c = (unsigned char*)str;
//   for (unsigned int i = 0 ; i < 512 ; ++i){
//     if (str[i] == '\0') break;
//     PushByteInBase64(c[i]);
//   }
// }
// /* -------------------------------------------------------------------------- */


// inline void Base64Writer::PushDoubleInBase64(double d){
//   //  DUMP("pushing double " << d << " as " << sizeof(double) << " bytes",DBG_ALL);

//   unsigned char * c = (unsigned char*)&d;
//   for (unsigned int i = 0 ; i < sizeof(double) ; ++i){
//     PushByteInBase64(c[i]);
//   }
// }
// /* -------------------------------------------------------------------------- */


inline void Base64Writer::PushByteInBase64(unsigned char c){
  //initialise les blocs
  //  DUMP("pushing byte " << (int) c << " at position " << n,DBG_ALL);

  if (n == 0){
    igroup[0]= 0;
    igroup[1]= 0;
    igroup[2]= 0;
  }
  igroup[n]= c;
  ++n;

  if(n == 3){
    dumpToBuffer();
  }
  nbBytes += 1;
}
/* -------------------------------------------------------------------------- */


inline void Base64Writer::finish(){
  if (n == 0) {
    return;
  }

  dumpToBuffer();
  linelength = 0;
}

/* -------------------------------------------------------------------------- */


inline void Base64Writer::dumpToBuffer(){
  if(n<3){
    igroup[2]= 0;
    if(n<2){
      igroup[1]= 0;
    }
  }

  //DUMP("premiere partie en base 64 : " << (igroup[0]>>2),DBG_ALL);
  int index = igroup[0]>>2;
  ogroup[0]= etable[index];
  //DUMP("deuxieme partie en base 64 : " << (((igroup[0]&3)<<4)|(igroup[1]>>4)),DBG_ALL);
  index = ((igroup[0]&3)<<4)|(igroup[1]>>4);
  ogroup[1]= etable[index];
  //DUMP("troisieme partie en base 64 : " << (((igroup[1]&0xF)<<2)|(igroup[2]>>6)),DBG_ALL);
  index = ((igroup[1]&0xF)<<2)|(igroup[2]>>6);
  ogroup[2]= etable[index];
  //DUMP("last partie en base 64 : " << (igroup[2]&0x3F),DBG_ALL);
  index = igroup[2]&0x3F;
  ogroup[3]= etable[index];

  if(n<3){
    ogroup[3]= '=';
    if(n<2){
      ogroup[2]= '=';
    }
  }

  for(int i= 0;i<4;i++){
    //DUMP("dumped to buffer " << ogroup[i],DBG_ALL);
    if (start == -1) {
      buffer.push_back(ogroup[i]);
    } else {
      buffer[start] = ogroup[i];
      ++start;
    }
  }

  //remise a zero du compteur
    n = 0;
}

/* inline void Base64Writer::ochar(int c) */
/* { */
/*   if(file.dumpchar(c)==-1){ */
/*       FATAL("error while writing to file (in compressed mode) ! no more space ?"); */
/*   } */
/*   linelength++; */
/* } */

/* -------------------------------------------------------------------------- */


inline int Base64Writer::Decode(char c0,char c1,char c2,char c3,
				 char * r0,char * r1,char * r2){
  auto d0 = dtable[0+c0];
  auto d1 = dtable[0+c1];
  auto d2 = dtable[0+c2];
  auto d3 = dtable[0+c3];

  //DUMP("d0 " << (int)d0 << " d1 " << (int)d1 << " d2 " << (int)d2 << " d3 " << (int)d3,DBG_ALL);

  // Decode the 3 bytes

  *r0 = (char)(((d0 << 2) & 0xFC) | ((d1 >> 4) & 0x03));
  *r1 = (char)(((d1 << 4) & 0xF0) | ((d2 >> 2) & 0x0F));
  *r2 = (char)(((d2 << 6) & 0xC0) | ((d3 >> 0) & 0x3F));

  //DUMP("r0 " << (int)*r0 << " r1 " << (int)*r1 << " r2 " << (int)*r2,DBG_ALL);

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
/* -------------------------------------------------------------------------- */



inline void Base64Writer::WriteHeader(){

  //  if (bflag == BASE64){
  //  char byteC[8];
  char byte[6];

  finish();

/*   long save_offset; */


/*   save_offset = file.tell(); */
/*   file.seek(header_offset,SEEK_SET); */
/*   //reread the 4 bytes precedently written */

/*   file.read(byteC,sizeof(char),4); */

/*   DUMP("la chaine saisie " << byteC[0] << " " << byteC[1] << " " << byteC[2] << " " << byteC[3]); */
/*   b64.Decode(byteC[0],byteC[0],byteC[0],byteC[0], */
/* 	     byte,byte+1,byte+2); */
  //DUMP("la chaine saisie " << buffer[0] << " " << buffer[1] << " " << buffer[2] << " " << buffer[3],DBG_ALL);

  Decode(buffer[0], buffer[0], buffer[0], buffer[0],
	 byte, byte + 1, byte + 2);

/*   file.read(byteC+4,sizeof(char),4); */
/*   DUMP("la chaine saisie " << byteC[4] << " " << byteC[5] << " " << byteC[6] << " " << byteC[7]); */
/*   int nb = b64.Decode(byteC[4],byteC[5],byteC[6],byteC[7], */
/* 		      byte+3,byte+4,byte+5); */

  //DUMP("la chaine saisie " << buffer[4] << " " << buffer[5] << " " << buffer[6] << " " << buffer[7],DBG_ALL);
  int nb = Decode(buffer[4],buffer[5],buffer[6],buffer[7],
		  byte + 3,byte + 4,byte + 5);


/*   //je me replace au debut du header */
/*   file.seek(header_offset,SEEK_SET); */
/*   //je viens de relire 6 octets */
/*   //les quatres premiers seulement sont a changer */


  //DUMP("placing number of writen bytes : " << nbBytes << " " << buffer.size(),DBG_ALL);

  auto temp = static_cast<int>(nbBytes);
  start = 0;
  push(temp);
  if (nb > 1) {
    push(byte[4]);
  }
  if (nb > 2) {
    push(byte[5]);
  }

  start = -1;
  nbBytes = temp;

  finish();
  DumpToFile();
}
/* -------------------------------------------------------------------------- */


inline void Base64Writer::DumpToFile(){
  file.write((char*)&buffer[0],buffer.size());
}
/* -------------------------------------------------------------------------- */


inline void Base64Writer::ClearBuffer(){
  buffer.clear();
  nbBytes = 0;
}
/* -------------------------------------------------------------------------- */


inline void Base64Writer::CreateHeader(){
  ClearBuffer();
  push<int>(0);
}

/* -------------------------------------------------------------------------- */

inline Base64Writer::Base64Writer(File & f):
  file(f){
  linelength = 0;
  InitBase64Stuff();
  n = 0;
  start = -1;
  igroup[0] = 0;
  igroup[1] = 0;
  igroup[2] = 0;
  ogroup[0] = 0;
  ogroup[1] = 0;
  ogroup[2] = 0;
  ogroup[3] = 0;
}
/* -------------------------------------------------------------------------- */

#if defined(__INTEL_COMPILER)
/// remark #981: operands are evaluated in unspecified order
#pragma warning ( pop )
#endif //defined(__INTEL_COMPILER)




}



#endif /* IOHELPER_BASE64_H_ */
