// -*- C++ -*-
// $Id: cppclass.h,v 1.1 2004/03/03 21:41:14 nrl Exp $
// 
// \file cppclass.h
// 
// Last modified: $Date: 2004/03/03 21:41:14 $
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
  CppClass();
  ~CppClass();
  
  CppClassMember const&  baseClass(int i) const { return baseList_[i]; }
  CppClassMember const&  member(int i) const { return memberList_[i]; }
  
  CppClass               instantiate(CppClassTemplateInstance const&) const;
  void                   copy(CppClass const&);
  CppClass&              operator=(CppClass const& t) { copy(t); return *this; }
  
private:
  std::vector<CppClassMember> baseList_;
  std::vector<CppClassMember> memberList_;
  
  void clear_();
};


class CppClassTemplateInstance : public CppTemplateInstance {
public:
  CppClassTemplateInstance() {}
  ~CppClassTemplateInstance() {}
  
  void    readFromHeaderSpec(std::string const& str) {}
};

#endif

