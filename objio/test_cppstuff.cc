// -*- C++ -*-
// $Id: test_cppstuff.cc,v 1.3 2004/03/08 13:20:44 guy Exp $
// 
// \file test_cppstuff.cc
// 

#include <iostream>
#include <vector>
#include <string>
#include "cpptype.h"
#include "cppclass.h"
#include "cppswigclassreader.h"

using namespace std;

int main()
{
  CppClass cl("class", "toto<T,U>", 2);
  CppClassMember m1("double", "toto");
  CppClassMember m2("double", "titi");
  CppClassMember m3("vector<string>", "titi");
  CppClassMember m4("vector<string>[12]", "titi");
  CppClassMember m5("list<vector<T> >", "glop");
  CppClassMember m6("T", "tshah_");
  CppClassMember m7("U", "tshuh_");
  
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
  
  CppSwigClassReader reader("sw.xml");
  int i;
  for(i=0;i<7;i++) {
    CppClass c1 = reader.read();
    c1.print(1);
    if(i==6) {
      vector<string> v;
      v.push_back("T"); v.push_back("U");
      c1.defineTemplateArgs(v);
      CppClass cc = c1.instantiate("BB<int,double>");
      cc.print(1);
    }
  }
  
  
}

