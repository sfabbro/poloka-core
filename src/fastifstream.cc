#include "fastifstream.h"

#include "cstdio"
//#include <string> // for isspace

/* about compression:
   we may compile this code so that it transparently reads gzipped files (or not).
   The burden is that with the zlib, the penalty is rather high:
   a 4.2 Mb (starlist) file is decoded in 
         - 0.4 s if read without zlib, 
	 - in 0.8 s if uncompressed and read through zlib,
         - in 1.8 s if compressed and read through zlib.

   So no that I have done the test, I think we should not use zlib.
   
*/

//#define USE_ZLIB


#ifdef USE_ZLIB
#include <zlib.h>
#endif


#include <iostream>

/* see fastifstream.h for the reason that lead to the implementation of this
   class.... and the considerable limitations of this implementation. */


#define LINE_BUFF_SIZE 8192

fastifstream::fastifstream(const char *name,const char *mode) 
{ 
#ifdef USE_ZLIB
  f = (FILE*) gzopen(name,mode); 
#else
  f = fopen(name,mode); 
#endif
  badbit = false; 
  line_buff = new char[LINE_BUFF_SIZE];
  if (f) get_next_line();
  else 
    {
      badbit = true;
      eofbit = true;
      return;
    }
  failbit = false;
  eofbit = file_eof();
}


bool fastifstream::file_eof() const
{
#ifdef USE_ZLIB
  return (gzeof((gzFile) f));
#else
  return (feof(f));
#endif
}

void fastifstream::get_next_line()
{
  if (file_eof())
    {
      eofbit = true;
      line_buff[0] = '\0';
      read_pointer = line_buff;
      return;
    }

#ifdef USE_ZLIB
  if (gzgets((gzFile)f, line_buff, LINE_BUFF_SIZE))
#else
  if (fgets(line_buff,LINE_BUFF_SIZE , f))
#endif
    {
      unsigned len = strlen(line_buff);
      if (len > LINE_BUFF_SIZE - 2)
	{
	  std::cout << "ERROR : (disaster) uncomplete line read in fastifstream.cc " << std::endl;
	  std::cerr << "ERROR : (disaster) uncomplete line read in fastifstream.cc " << std::endl;
	  abort();
	}
      if (len>=1 && line_buff[len-1] == '\n') line_buff[len-1] = '\0'; 
      read_pointer = line_buff;
    }
  else badbit = true; // not sure about the meaning of badbit in ifstream
}

void fastifstream::skip_spaces() // in ifstream, before extraction of (e.g.) a float
{
  while (isspace(*read_pointer)) read_pointer++;
  if (*read_pointer == '\0' && !file_eof()) 
    {
      get_next_line();
      skip_spaces();
    }
}


void fastifstream::close()
{
#ifdef USE_ZLIB
  if (f) gzclose((gzFile) f);
#else
  if (f) fclose(f);
#endif
  f = NULL;
}


fastifstream::~fastifstream()
{ 
  close();
  delete [] line_buff;
}  // fstream objects close the file in destructorsconst char *mode = "r"
