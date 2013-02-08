// -*- C++ -*-
// 
// A few helper functions to read/write class members
// and STL containers to an ASCII stream.
// 
#ifndef IOHELPERS_H
#define IOHELPERS_H


#include <iostream>
#include <string>



template<class T>
void write(std::ostream& os, std::vector<T> const& vec)
{
  size_t i;
  os << " " << vec.size() << endl;
  for(i=0;i<vec.size();i++)
    os << " - " << vec[i] << endl;
}


template<class T>
void read(std::istream& is, std::vector<T>& vec)
{
  vec.clear();
  
  size_t i, sz;
  T value;
  std::string tmp_str;
  is >> sz;
  for(i=0;i<sz;i++)
    {
      is >> tmp_str >> value;
      vec.push_back(value);
    }
}


template<class T>
std::ostream& operator<<(std::ostream& os, std::vector<T> const& vec)
{
  write(os, vec);
  return os;
}


template<class T>
std::istream& operator>>(std::istream& is, std::vector<T>& vec)
{
  read(is, vec);
  return is;
}


template<class T>
void write_member(std::ostream& os, 
		  std::string const& name, T const& value)
{
  os << std::endl << name << " " << value;
}


template<class T>
void read_member(std::istream& is, string& name, T& val)
{
  is >> name;
  is >> val;
}


#endif


