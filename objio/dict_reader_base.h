// -*- C++ -*-
// 
// file dict_reader_base.h
// 
// 
#ifndef DICT_READER_BASE_H
#define DICT_READER_BASE_H

#include <string>

#include "dict.h"


class dict_reader_base {
public:
  dict_reader_base(std::string const& filename) : filename_(filename) {}
  ~dict_reader_base() {}
  
  virtual unsigned int size() const=0;
  std::string          filename() const { return filename_; }
  inline virtual void  read(dict& d) const;
  
protected:
  std::string filename_;
  
  virtual void next_() const=0;
  virtual void read_header_(std::string&, unsigned int&, std::string&,
			    bool&, bool&) const=0;
  virtual void read_baselist_(std::vector<std::string>&) const=0;
  virtual void read_memberList_(std::vector<dict::member_t>&) const=0;
};




////////////////////////////// INLINED STUFF //////////////////////////////

void dict_reader_base::read(dict& d) const
{
  read_header_(d.name_, d.version_, d.kind_, 
	       d.isTemplate_, d.isPersistent_);
  d.realName_=d.symName_=d.name_;
  read_baselist_(d.baseList_);
  read_memberList_(d.memberList_);
  next_();
}


#endif

