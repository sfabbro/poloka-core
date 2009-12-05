#ifndef FASTIFSTREAM__H
#define FASTIFSTREAM__H



//#define USE_STL_IFSTREAM
#ifndef USE_STL_IFSTREAM

#include <cstdio>
#include <string>
#include <cctype>
#include <string.h> // for strcpy
#include <stdlib.h> // for strtol strtoul

/*! This is a Mickey-Mouse implementation of ifstream over core C I/O's.
The reason is that it is faster that actual ifstream: by a factor of 2
for gcc version 4.1 and a factor of 5 (!) for gcc 3.2 . However, this implementation
is not only incomplete (it misses a lot of formats), it is also wrong!
In particular the "char" reading and unget routines work for what StarList::read does,
but not necessarily in general. For example, extremely long lines are not handled.
*/


class fastifstream {
private :
  FILE *f;
  char last_char_read;
  char *line_buff;
  char *read_pointer;
  bool eofbit;
  bool badbit;
  bool failbit;

private :
  void get_next_line();

  bool file_eof() const;

  void skip_spaces(); // in ifstream, before extraction of (e.g.) a float

public :
  fastifstream(const char *name, const char *mode = "r");

  operator bool() const { return (f && !eofbit && !badbit);}
  

  bool eof() const {return eofbit;}
  bool bad() const { return badbit;}
  bool fail() const {return failbit;}
  //! in istream, good is not the complement of bad....
  bool good() const {return f && !badbit && !failbit && !eofbit;}

  int rdstate() const { return eofbit*4+2*failbit+badbit;}

  fastifstream& operator >> (double &d) 
  { 
    skip_spaces();
    char *end=NULL;
    d = strtod(read_pointer,&end); 
    if (end==read_pointer)
      { 
	failbit=true;
	return (*this);
      }
    else {read_pointer = end;}
    return *this;
  }


  fastifstream& operator >> (int &i) 
  { 
    skip_spaces();
    char *end=NULL;
    i = strtol(read_pointer,&end,10); 
    if (end==read_pointer)
      { 
	failbit=true;
	return (*this);
      }
    else {read_pointer = end;}
    return *this;
  }

  fastifstream& operator >> (unsigned &i) 
  { 
    skip_spaces();
    char *end=NULL;
    i = strtoul(read_pointer,&end,10); 
    if (end==read_pointer)
      { 
	failbit=true;
	return (*this);
      }
    else {read_pointer = end;}
    return *this;
  }


  fastifstream& operator >>   (char &c) 
  { 
    skip_spaces();
    c = *read_pointer; 
    read_pointer++;
    return *this;
  }
  
  fastifstream& getline(char *buff, size_t len)
  { 
    strncpy(buff,read_pointer,len); 
    get_next_line();
    return *this;
  }

  fastifstream& operator >>   (std::string &s)
  {
    skip_spaces();
    s.clear();
    while (*read_pointer != '\0' && !isspace(*read_pointer)) 
      {
	s.push_back(*read_pointer);
	read_pointer++;
      }
    return *this;
  }
    
    


  void close();

  void unget() {read_pointer--;}

  ~fastifstream();

 private:
  // forbid copies
  fastifstream(const fastifstream&);
  fastifstream& operator = (const fastifstream &);

};

#else

#include <fstream>

class fastifstream : public std::ifstream
{
 public :
  fastifstream(const std::string &FileName) : std::ifstream(FileName.c_str()) {};
};
#endif





#endif /* FASTIFSTREAM__H */
