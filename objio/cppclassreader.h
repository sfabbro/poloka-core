// -*- C++ -*-
// $Id: cppclassreader.h,v 1.1 2004/03/08 11:50:41 nrl Exp $
// \file cppclassreader.h
// 
// 
#ifndef CPPCLASSREADER_H
#define CPPCLASSREADER_H

#include <string>
#include "cppclass.h"


class CppClassReader {
public:
  CppClassReader() {}
  CppClassReader(std::string const& filename) : filename_(filename) {}
  virtual ~CppClassReader() {}
  
  virtual int          size() const=0;
  virtual CppClass     read() const=0;
  
protected:
  std::string filename_;
};


#endif

