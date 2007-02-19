#ifndef FASTIFSTREAM__H
#define FASTIFSTREAM__H



//#define USE_STL_IFSTREAM
#ifndef USE_STL_IFSTREAM

#include <cstdio>
#include <string>

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
  char line_buff[8192];
  char *read_pointer;
  bool badbit;

private :
  void get_next_line()
  {
    fgets(line_buff, 8192, f);
    unsigned len = strlen(line_buff);
    if (len>8190)
      {
	cout << " disaster : uncomplete line read int fastifstream.h " << endl;
	abort();
      }
    if (len>=1 && line_buff[len-1] == '\n') line_buff[len-1] = '\0'; 
    read_pointer = line_buff;
  }

  void skip_spaces() // in ifstream, before extraction of (e.g.) a float
  {
    while (isspace(*read_pointer)) read_pointer++;
    if (*read_pointer == '\0' && !feof(f)) 
      {
	get_next_line();
	skip_spaces();
      }
  }

public :
  fastifstream(const char *name) 
    { f = fopen(name,"r"); badbit = false;
      if (f) get_next_line();}

  operator bool() const { return (f && !feof(f) && !badbit);}
  
  bool bad() const { return badbit;}

  fastifstream& operator >> (double &d) 
  { 
    skip_spaces();
    char *end=NULL;
    d = strtod(read_pointer,&end); 
    if (end==read_pointer)
      { 
	badbit=true;
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
	badbit=true;
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
	badbit=true;
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
    while (!isspace(*read_pointer)) 
      {
	s.push_back(*read_pointer);
	read_pointer++;
      }
    return *this;
  }
    
    


  void close() { fclose(f); f = NULL;} 

  void unget() {read_pointer--;}

  ~fastifstream() { fclose(f);}  // fstream objects close the file in destructors

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
