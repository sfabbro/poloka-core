// -*- C++ -*-
// $Id: cppclass.cc,v 1.4 2004/03/08 13:20:01 guy Exp $
// \file cppclass.cc
// 
// 
#include <iostream>
#include <string>
#include <algorithm>
#include "cppclass.h"

using namespace std;


CppClass::CppClass(std::string const& kind, std::string const& classname, int version,
		   bool isTemplate)
  : CppType(classname), version_(version), kind_(kind)
{
  isTemplate_=isTemplate;
}


CppClass::~CppClass()
{
}


CppClass CppClass::instantiate(std::string const& str) const
{
  CppClass ret(kind_,cppTypeName_,version_);
  if( !isTemplate() ) {
    cout << "CppClass::instantiate: class is not template !"
	 << endl;
    ret = *this;
    return ret;
  }
  
  // build the CppTemplateInstance
  CppTemplateInstance ti = readFromHeaderSpec(str);
  if(ti.size()==0) {
    ret = *this;
    return ret;
  }
  
  // instantiate the type itself
  CppType t = CppType::instantiate(ti);
  ret.CppType::copy(t);
  
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
  
  CppType::copy((CppType const&)cp);
  version_ = cp.version_;
  kind_ = cp.kind_;
  std::copy(cp.baseList_.begin(),cp.baseList_.end(),
	    back_inserter(baseList_));
  std::copy(cp.memberList_.begin(),cp.memberList_.end(),
	    back_inserter(memberList_));
}


void CppClass::clear_()
{
  version_=0;
  kind_="";
  baseList_.clear();
  memberList_.clear();
}



CppTemplateInstance CppClass::readFromHeaderSpec(string const& str) const
{
  CppTemplateInstance ret;
  
  //  vector<string> tok;
  //  CppType::tokenize(str,tok,"<> \t,",true);
  CppType t(str);
  
  if(t.nTemplateArgs()!=nTemplateArgs()) {
    cout << "CppClass::readFromHeaderSpec(" << str << "): ERROR "
	 << "template argument mismatch:"  << endl
	 << " class: " << nTemplateArgs() << " template args != "
	 << t.nTemplateArgs()
	 << endl;
    return ret;
  }
  

  // now, build the CppTemplateInstance
  int i,sz=t.nTemplateArgs();
  for(i=0;i<sz;i++)
    ret.addInstance(templateArg(i),t.templateArg(i));
  
  return ret;
}


void CppClass::addMember(CppClassMember const& m)
{
  memberList_.push_back(m);
}


void CppClass::addBaseClass(CppClassMember const& b)
{
  baseList_.push_back(b);
}


void CppClass::defineTemplateArgs(std::vector<std::string> const& args)
{
  if(args.size()==0) return;
  if(!isTemplate_) {
    cout << "CppClass::defineTemplateArgs() class is not template!" 
	 << endl;
    return;
  }
  templateArgs_.clear();
  std::copy(args.begin(), args.end(),
	    std::back_inserter(templateArgs_));
  
  int i,sz=args.size();
  if(sz>0) {
    tok_.push_back("<");
    baseTok_.push_back("<");
    tok_.push_back(args[0]);
    baseTok_.push_back(args[0]);
  }
  for(i=1;i<sz;i++) {
    tok_.push_back(",");
    baseTok_.push_back(",");
    tok_.push_back(args[i]);
    baseTok_.push_back(args[i]);
  }
  if(sz>0) {
    tok_.push_back(">");
    baseTok_.push_back(">");
  }
  
  cppTypeName_ = buildTypeString(tok_);
  symbolicCppTypeName_ = buildTypeString(baseTok_);
}


void CppClass::print(int verbosity) const
{
  if( isTemplate() && !wasInstantiated_ )
    std::cout << "T ";
  else
    std::cout << "  ";
  
  std::cout.flags(ios::left); std::cout.width(10); 
  std::cout << kind_.c_str();
  std::cout << " " << cppTypeName_ << "[" 
	    << version_ << "] ";
  
  int i,sz=baseList_.size();
  if(sz!=0) {
    std::cout << "{ ";
    for(i=0;i<sz;i++)
      std::cout << baseList_[i].type().cppTypeName() << " ";
    std::cout << "}";
  }
  
  if( isTemplate() ) {
    std::cout << "  < ";
    sz=templateArgs_.size();
    for(i=0;i<sz;i++)
      std::cout << templateArgs_[i] << " ";
    std::cout << ">";
  }
  std::cout << std::endl;  
  
  if(verbosity>0) {
    sz=memberList_.size();
    for(i=0;i<sz;i++) 
      memberList_[i].print();
  }
}







  // first, check the class name...
  //  string cname = cppTypeName_.substr(0,cppTypeName_.find_first_of(" <\t"));
  //  if(cname != tok[0]) {
  //    cout << "CppClass::readFromHeaderSpec("
  //	 << str << ") required class name "
  //	 << tok[0]
  //	 << " does not match the CppClass Name "
  //	 << cppTypeName_ << endl;
  //  }
  
  //  // parse the tokens to find the template args B<int,map<string,string> >
  //  int i,depth=0,sz=tok.size();
  //  string tpl_arg;
  //  vector<string> tpl_arg_vec;
  //  for(i=0;i<sz;i++) {
  //    
  //    if(tok[i]=="<") {
  //      if(depth==0) {
  //	if(tpl_arg.size()>0)
  //	  tpl_arg_vec.push_back(tpl_arg);
  //	tpl_arg="";
  //      }
  //      if(depth>0)
  //	tpl_arg = tpl_arg + tok[i];
  //      depth++;
  //      continue;
  //    }
  //    
  //    if(tok[i]==">") {
  //      if(depth==1) {
  //	tpl_arg_vec.push_back(tpl_arg);
  //	tpl_arg="";
  //      }
  //      if(depth>1)
  //	tpl_arg = tpl_arg + tok[i];
  //      depth--;
  //      continue;
  //    }
  //    
  //    if(tok[i]==",") {
  //      if(depth==1) {
  //	tpl_arg_vec.push_back(tpl_arg);
  //	tpl_arg="";
  //      }
  //      if(depth>1) 
  //	tpl_arg = tpl_arg + tok[i];
  //      continue;
  //    }
  //    
  //    if(depth>=1)
  //      tpl_arg = tpl_arg + tok[i];
  //  }
  
  //  sz=templateArgs_.size();
  //  for(i=0;i<sz;i++) {
  //    string::size_type start = tpl_arg_vec[i].find_first_not_of(" \t");
  //    string::size_type stop  = tpl_arg_vec[i].find_last_not_of(" \t");
  //    tpl_arg_vec[i] = tpl_arg_vec[i].substr(start,stop-start);
  //    cout << tpl_arg_vec[i] << endl;
  //  }
  
  //  if(tpl_arg_vec.size() != templateArgs_.size()) {
  //    cout << "warning ! something nasty going on!" << endl;
  //    return ret;
  //  }
  
  // now, build the CppTemplateInstance
  //  for(i=0;i<sz;i++)
  //    ret.addInstance(templateArgs_[i],tpl_arg_vec[i]);
  
