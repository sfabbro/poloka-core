// -*- C++ -*-
// $Id: cppclass.cc,v 1.1 2004/03/03 21:41:14 nrl Exp $
// \file cppclass.cc
// 
// 
#include <algorithm>
#include "cppclass.h"



CppClass::CppClass()
  : CppType()
{
}


CppClass::~CppClass()
{
}


CppClass CppClass::instantiate(CppClassTemplateInstance const& ti) const
{
  CppClass ret;
  ret=*this;
  
  int i,sz=baseList_.size();
  for(i=0;i<sz;i++) {
    CppClassMember m;
    m.name() = baseList_[i].name();
    m.type() = baseList_[i].type().instantiate(ti);
    ret.baseList_.push_back(m);
  }
  sz=memberList_.size();
  for(i=0;i<sz;i++) {
    CppClassMember m;
    m.name() = memberList_[i].name();
    m.type() = memberList_[i].type().instantiate(ti);
    ret.memberList_.push_back(m);
  }
  
  return ret;
}


void CppClass::copy(CppClass const& cp)
{
  clear_();
  
  std::copy(cp.baseList_.begin(),cp.baseList_.end(),
	    back_inserter(baseList_));
  std::copy(cp.memberList_.begin(),cp.memberList_.end(),
	    back_inserter(memberList_));
}


void CppClass::clear_()
{
  baseList_.clear();
  memberList_.clear();
}
