// -*- C++ -*-
// $Id: cppclassmember.cc,v 1.3 2004/03/06 23:15:42 nrl Exp $
// 
// \file cppclassmember.cc
// 
#include "cppclassmember.h"


CppClassMember::CppClassMember()
  : name_("")
{
}


CppClassMember::CppClassMember(CppType const& decl, std::string const& name)
  : name_(name), type_(decl)
{
}


CppClassMember::CppClassMember(const std::string& decl, const std::string& name)
  : name_(name), type_(decl)
{
}


void CppClassMember::copy(CppClassMember const& cm)
{
  clear();
  name_=cm.name_;
  type_=cm.type_;
}


void CppClassMember::print() const
{
  if( type_.isPersistent() )
    std::cout << "   + ";
  else
    std::cout << "   - ";
  std::cout.flags(ios::left); std::cout.width(30); 
  std::cout << type_.cppTypeName().c_str();
  
  std::cout.flags(ios::left);// std::cout.width(45); 
  std::cout << name_.c_str();
  //  int i,sz=arraySize.size();
  //  for(i=0;i<sz;i++)
  //    std::cout << "[" << arraySize[i] << "]";
  
  std::cout << std::endl;
}


void CppClassMember::clear()
{
  name_="";
  type_.clear();
}
