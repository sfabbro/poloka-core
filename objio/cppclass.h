// -*- C++ -*-
// $Id: cppclass.h,v 1.2 2004/03/04 17:49:06 nrl Exp $
// 
// \file cppclass.h
// 
// Last modified: $Date: 2004/03/04 17:49:06 $
// By:            $Author: nrl $
// 
#ifndef CPPCLASS_H
#define CPPCLASS_H

#include <vector>

#include "cpptype.h"
#include "cppclassmember.h"


class CppClassTemplateInstance;


class CppClass : public CppType {
public:
  CppClass(std::string const& kind, std::string const& classname, int version);
  ~CppClass();
  
  CppClassMember const&  baseClass(int i) const { return baseList_[i]; }
  CppClassMember const&  member(int i) const { return memberList_[i]; }
  
  int                    version() const { return version_; } 
  CppClass               instantiate(std::string const&) const;
  
  void                   copy(CppClass const&);
  CppClass&              operator=(CppClass const& t) { copy(t); return *this; }

  //  void                   copy(CppType const&);
  //  CppClass&              operator=(CppType const& t) { copy(t); return *this; }
  
  CppTemplateInstance    readFromHeaderSpec(std::string const& str) const;
  void                   addMember(CppClassMember const&);
  void                   addBaseClass(CppClassMember const&);
  
  void                   print(int verbosity=0) const;
  
private:
  int         version_;
  std::string kind_;
  std::vector<CppClassMember> baseList_;
  std::vector<CppClassMember> memberList_;
  
  void clear_();
};


#endif

