// -*- C++ -*-
// $Id: test_cppstuff.cc,v 1.1 2004/03/04 17:50:04 nrl Exp $
// 
// \file test_cppstuff.cc
// 
#include "cpptype.h"
#include "cppclass.h"



int main()
{
  CppClass cl("class", "toto<T,U>", 2);
  CppClassMember m1(CppType("double"), "toto");
  CppClassMember m2(CppType("double"), "titi");
  CppClassMember m3(CppType("vector<string>"), "titi");
  CppClassMember m4(CppType("vector<string>[12]"), "titi");
  CppClassMember m5(CppType("list<vector<T> >"), "glop");
  CppClassMember m6(CppType("T"), "tshah_");
  CppClassMember m7(CppType("U"), "tshuh_");
  
  cl.addMember(m1);
  cl.addMember(m2);
  cl.addMember(m3);
  cl.addMember(m4);
  cl.addBaseClass(m5);
  cl.addMember(m6);
  cl.addMember(m7);
  cl.print(1);
  
  CppClass cl2 = cl.instantiate("toto<vector<map<int,double> >, int >");
  cl2.print(1);
}

