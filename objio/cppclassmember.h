// -*- C++ -*-
// $Id: cppclassmember.h,v 1.1 2004/03/03 21:41:15 nrl Exp $
// 
// \file cppclassmember.h
// 
// Last modified: $Date: 2004/03/03 21:41:15 $
// By:            $Author: nrl $
// 
#ifndef CPPCLASSMEMBER_H
#define CPPCLASSMEMBER_H

#include <string>

#include "cpptype.h"


class CppClassMember {
public:
  CppClassMember();
  ~CppClassMember() {}
  
  std::string const& name() const { return name_; }
  CppType const&     type() const { return type_; }
  std::string&       name()       { return name_; }
  CppType&           type()       { return type_; }
  
  void               copy(CppClassMember const&);
  CppClassMember&    operator=(CppClassMember const& m) { copy(m); return *this; }
  
private:
  std::string name_;
  CppType     type_;
};


#endif

