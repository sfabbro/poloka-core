#ifndef FASTIFSTREAM__H
#define FASTIFSTREAM__H



//#define USE_STL_IFSTREAM
#ifndef USE_STL_IFSTREAM

#include <cstdio>

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
    if (len>1 && line_buff[len-1] == '\n') line_buff[len-1] = '\0'; 
    read_pointer = line_buff;
  }

  void skip_spaces() // in ifstream, after a read of (e.g.) a float
  {
    while (*read_pointer == ' ') read_pointer++;
  }

public :
  fastifstream(const char *name) 
  { f = fopen(name,"r"); if (f) get_next_line();  }

  operator bool() const { return (f && !feof(f));}

  fastifstream& operator >> (double &d) 
  { 
    char *end=NULL;
    d = strtod(read_pointer,&end); 
    if (end==read_pointer) {get_next_line(); return ( (*this) >> d);}
    else {read_pointer = end;}
    skip_spaces();
    return *this;
  }
  fastifstream& operator >> (float &d) 
  { 
    char *end=NULL;
    d = strtod(read_pointer,&end); 
    if (end==read_pointer) {get_next_line(); return ( (*this) >> d);}
    else {read_pointer = end;}
    skip_spaces();
    return *this;
  }

  fastifstream& operator >> (int &i) 
  { 
    char *end=NULL;
    i = strtol(read_pointer,&end,10); 
    if (end==read_pointer) {get_next_line(); return ( (*this) >> i);}
    else {read_pointer = end;}
    skip_spaces();
    return *this;
  }

  fastifstream& operator >> (unsigned &i) 
  { 
    char *end=NULL;
    i = strtoul(read_pointer,&end,10); 
    if (end==read_pointer) {get_next_line(); return ( (*this) >> i);}
    else {read_pointer = end;}
    skip_spaces();
    return *this;
  }


  fastifstream& operator >>   (char &c) 
  { 
    if (*read_pointer == '\0') get_next_line();
    c = *read_pointer; 
    read_pointer++;
    return *this;
  }
#include <string.h>
 fastifstream& operator >> (std::string& s)
  {
    while (*read_pointer == ' ') read_pointer++;
    s.clear();
    do
      {
        s.push_back(*read_pointer);
        read_pointer++;
      }
    while (*read_pointer != ' ' && *read_pointer != '\0');
    if (*read_pointer == '\0') get_next_line();
    return *this;
  }
  
  fastifstream& getline(char *buff, size_t len)
  { 
    strncpy(buff,read_pointer,len); 
    get_next_line();
    return *this;
  }

  void close() { fclose(f); f = NULL;}

  void unget() {read_pointer--;}

  ~fastifstream() { fclose(f);}

};

#else

#include <fstream>

class fastifstream : public ifstream
{
 public :
  fastifstream(const string &FileName) : ifstream(FileName.c_str()) {};
};
#endif





#endif /* FASTIFSTREAM__H */
