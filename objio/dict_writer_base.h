// -*- C++ -*-
// 
// file dict_writer_base.h
// 
// 
#ifndef DICT_WRITER_BASE_H
#define DICT_WRITER_BASE_H


class dict; 

class dict_writer_base {
public:
  dict_writer_base(std::string const& filename) : filename_(filename) {}
  ~dict_writer_base() {}
  
  std::string       filename() const { return filename_; }
  virtual void      write(dict const&)=0;

protected:
  std::string filename_;
};


#endif


