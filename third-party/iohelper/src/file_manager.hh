/**
 * @file   file_manager.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Thu Dec 06 2012
 *
 * @brief  file manager header
 *
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef IOHELPER_FILE_MANAGER_H_
#define IOHELPER_FILE_MANAGER_H_
/* -------------------------------------------------------------------------- */
#include <zlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include "iohelper_common.hh"
/* -------------------------------------------------------------------------- */

namespace iohelper {

template <class charT, class Traits=std::char_traits<charT> >
class GZfstream : public std::basic_fstream<charT,Traits> {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:

  GZfstream();
  GZfstream(const std::string & fname,
		     std::fstream::openmode mode = std::fstream::out,
		     bool compr=false);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  ///! opening methods
  inline void open(const std::string & name,
		   std::fstream::openmode mode = std::fstream::out,
		   bool compr= false);
  // inline void close();

  // // //! writing methods
  // template <typename T> GZfstream & operator << (const T & v);
  // inline GZfstream & operator << (std::ostream& (*op)(std::ostream&)){};
  // inline GZfstream & write(const void * buffer,std::streamsize n);
  // inline GZfstream & flush();
  // // //! reading methods
  // template <typename T> GZfstream & operator >> (T & v);
  // inline GZfstream & read(void * buffer,std::streamsize n);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */


 private:

// #ifdef USING_ZLIB
//   gzFile gzfile;
// #endif

//   int compressed;

};


// /* -------------------------------------------------------------------------- */


template <class charT, class Traits>
 inline GZfstream<charT,Traits>::GZfstream():
   std::basic_fstream<charT,Traits>(){
//     compressed = 0;
// #ifdef USING_ZLIB
//     gzfile = NULL;
// #endif
}

// /* -------------------------------------------------------------------------- */

template <class charT, class Traits>
inline GZfstream<charT, Traits>::GZfstream(const std::string & fname,
                                           std::fstream::openmode mode,
                                           bool /*unused*/)
    : std::basic_fstream<charT, Traits>(fname.c_str(), mode) {
  //  compressed = compr;
// #ifdef USING_ZLIB
//   gzfile = NULL;
// #endif
//  open(fname,mode,compressed);
}

// /* -------------------------------------------------------------------------- */

// // inline void GZfstream::printf(const std::string & formated, ...){
// //   if (!opened) {
// //     FATAL("Warning : file not opened " << name);
// //   }


// //   //  if (!opened) return;

// //   if (_file == NULL) FATAL("fichier non ouvert mais ce n'est pas normal");

// //   char buf[512];

// //   va_list list;
// //   va_start(list,formated);
// //   int len = vsprintf(buf,formated.c_str(),list);
// //   va_end(list);

// // #ifdef USING_ZLIB
// //     if (compressed){
// //       gzwrite(gzfile,buf,len);
// //       return;
// //     }
// //   else
// // #endif
// //     fwrite(buf,len,1,_file);

// // }

// /* -------------------------------------------------------------------------- */


// // inline void GZfstream::gets(std::string & buf){
// //   char * ret;
// //   const UInt len = 255;
// //   char buffer[len] = "";

// //   if (!opened) {
// //     FATAL("Warning : file not opened " << name);
// //   }
// //   if (_file == NULL) FATAL("file not opened: exit");

// //   DUMP("read to buf (" << buf << ") at most " << len << " characters");
// // #ifdef USING_ZLIB
// //     if (compressed){
// //       ret = gzgets(gzfile,buffer,len);
// //     }
// //   else
// // #endif
// //     ret = fgets(buffer,len,_file);

// //     if (ret == NULL) throw;

// //     DUMP("read to buf (" << buf << ") at most " << len << " characters");

// //     buf = buffer;
// // }

// /* -------------------------------------------------------------------------- */


// template <class charT, class Traits>
// inline void GZfstream<charT,Traits>::close(){

// #ifdef USING_ZLIB
//   if (this->is_open() && compressed){
//     gzclose(gzfile);
//     gzfile = NULL;
//   }
// #endif
// }

// /* -------------------------------------------------------------------------- */

template <class charT, class Traits>
inline void GZfstream<charT, Traits>::open(const std::string & fname,
                                           std::fstream::openmode mode,
                                           bool /*unused*/) {
  std::basic_fstream<charT,Traits>::open(fname.c_str(),mode);

//   if (this->is_open()) this->close();

//   compressed = compr;
//   this->open(fname.c_str(),mode);
//   if (!this->is_open()) FATAL("Could not open file "<< fname);

// #ifdef USING_ZLIB
//   if (compressed){
//     std::stringstream _mode;
//     if (mode & std::fstream::in) _mode  << "r";
//     if (mode & std::fstream::out) _mode << "w";
//     if (mode & std::fstream::app) _mode << "a";
//     if (mode & std::fstream::binary) _mode << "b";

//     std::cerr << typeid(*this->rdbuf()).name() << std::endl;
//     int fd = this->rdbuf()->_M_file->fd();

//     //gzfile = gzdopen(this->rdbuf()->fd(),_mode.str().c_str());
//     throw;
//   }
// #endif

//   DUMP("file " << name << " opened , compressed = " << compressed);
}

// /* -------------------------------------------------------------------------- */


// // inline int GZfstream::dumpchar(int c){
// //   if (!opened) {
// //     FATAL("Warning : file not opened " << name);
// //   }

// //   //  if (!opened) return EOF;

// // #ifdef USING_ZLIB
// //   if (compressed){
// //     int res = gzputc(gzfile,c);
// //     //DUMP("file opened in compressed form " << name);
// //     if (c != res)
// //       FATAL("j'ai pas ecrit ce que je voulais (compressed) " << res);
// //     return res;
// //   }
// //   else
// // #endif
// //     if (c != putc(c,_file))
// //       FATAL("j'ai pas ecrit ce que je voulais");
// //   return c;
// // }
// /* -------------------------------------------------------------------------- */



// // inline int GZfstream::seek(int offset,int set){
// //   if (!opened) {
// //     FATAL("Warning : file not opened " << name);
// //   }


// // #ifdef USING_ZLIB
// //   if (compressed){
// //     return gzseek(gzfile, offset,set);
// //   }
// //   else
// // #endif
// //     return fseek(_file,offset,set);
// // }

// // /* -------------------------------------------------------------------------- */


// // inline int GZfstream::tell(){
// //   if (!opened) {
// //     FATAL("Warning : file not opened " << name);
// //   }

// // #ifdef USING_ZLIB
// //   if (compressed){
// //     return gztell(gzfile);
// //   }
// //   else
// // #endif

// //     return ftell(_file);
// // }

// // /* -------------------------------------------------------------------------- */

// //template <class charT, class Traits>
// //inline GZfstream<charT,Traits> & GZfstream<charT,Traits>::read(void * buffer,
// //							       std::streamsize size){
// //   if (!opened) {
// //     FATAL("Warning : file not opened " << name);
// //   }


// // #ifdef USING_ZLIB
// //   if (compressed){
// //     return gzread(gzfile,buffer,size*number);
// //   }
// //   else
// // #endif
// //     return fread(buffer,size,number,_file);
// //}

// /* -------------------------------------------------------------------------- */


// // template <class charT, class Traits>
// // inline GZfstream<charT,Traits> & GZfstream<charT,Traits>::write(const void * buffer,
// // 								std::streamsize size){
// //   //  if (!this->is_open()) throw std::ios_base::failure("file not opened");

// //   Int nwrite;
// // #ifdef USING_ZLIB
// //   if (compressed) nwrite = gzwrite(gzfile,buffer,size);
// //   else
// // #endif
// //     this->write(buffer,size);

// //   //  if (nwrite == 0) throw std::ios_base::failure("could not write any byte");

// //   return *this;
// // }

// /* -------------------------------------------------------------------------- */


// // template <class charT, class Traits>
// // inline GZfstream<charT,Traits> & GZfstream<charT,Traits>::flush(){
// //   //  if (!this->is_open()) throw std::ios_base::failure("file not opened");

// // #ifdef USING_ZLIB
// //   if (compressed){
// //      gzflush(gzfile,Z_SYNC_FLUSH);
// //   }
// //   else
// // #endif
// //     this->flush();
// // }

using File = GZfstream<char>;

}




#endif /* IOHELPER_FILE_MANAGER_H_ */
