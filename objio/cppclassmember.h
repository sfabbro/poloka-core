// -*- C++ -*-
// 
// \file cppclassmember.h
// 
// 
#ifndef CPPCLASSMEMBER_H
#define CPPCLASSMEMBER_H

#include <string>

#include "cpptype.h"


class CppClassMember {
public:
  CppClassMember();
  CppClassMember(CppType const& decl, std::string const& name);
  CppClassMember(const std::string& decl, const std::string& name);
  ~CppClassMember() {}
  
  std::string const& name() const { return name_; }
  CppType const&     type() const { return type_; }
  std::string&       name()       { return name_; }
  CppType&           type()       { return type_; }
  
  void               copy(CppClassMember const&);
  CppClassMember&    operator=(CppClassMember const& m) { copy(m); return *this; }
  
  void               print() const;
  
  void               clear();
private:
  std::string name_;
  CppType     type_;
};


#endif

