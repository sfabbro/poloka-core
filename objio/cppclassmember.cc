// -*- C++ -*-
// $Id: cppclassmember.cc,v 1.1 2004/03/03 21:41:15 nrl Exp $
// 
// \file cppclassmember.cc
// 
#include "cppclassmember.h"


CppClassMember::CppClassMember()
  : name_("")
{
}


void CppClassMember::copy(CppClassMember const& cm)
{
  name_=cm.name_;
  type_=cm.type_;
}

