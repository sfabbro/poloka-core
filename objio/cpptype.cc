// -*- C++ -*-
// $Id: cpptype.cc,v 1.1 2004/03/03 21:41:15 nrl Exp $
// \file cpptype.cc
// 
// 
#include <algorithm>
#include "cpptype.h"

using namespace std;


CppType::CppType()
  : isConst_(false),isStatic_(false),
    isAPointer_(false), // isASTLContainer_(false),
    isAComplexType_(false), isAReference_(false),
    wasInstantiated_(false)
{
}


CppType::~CppType()
{
}


bool CppType::isPersistent() const
{
  if(cppTypeName_=="") return false;
  if(isConst_) return false;
  if(isStatic_) return false;
  if(isAPointer_) return false;
  if(isAReference_) return false;
  return true;
}


CppType CppType::instantiate(CppTemplateInstance const& ti) const
{
  CppType ret;
  ret=*this;
  
  if(wasInstantiated_) {
    cout << "CppType::instantiate() ERROR"
	 << " class was already instantiated..." 
	 << endl;
    return ret;
  }
  
  if(!isATemplate()) {
    cout << "CppType::instantiate() ERROR"
	 << " class is not template!" 
	 << endl;
    return ret;
  }
  
  int i,idx,sz=tok_.size();
  for(i=0;i<sz;i++) {
    if(ti.hasSymName(tok_[i])) 
      ret.tok_[i]=ti.realName(tok_[i]);
  }
  
  sz=baseTok_.size();
  for(i=0;i<sz;i++) {
    if(ti.hasSymName(baseTok_[i])) 
      ret.baseTok_[i]=ti.realName(baseTok_[i]);
  }
  
  // now, we recompute the instantiated names
  ret.symbolicCppTypeName_=cppTypeName_;
  ret.symbolicBaseCppTypeName_=baseCppTypeName_;
  ret.cppTypeName_="";
  ret.baseCppTypeName_="";
  sz=ret.tok_.size();
  for(i=0;i<sz;i++) 
    ret.cppTypeName_ = ret.cppTypeName_ + ret.tok_[i] + " ";
  
  sz=ret.baseTok_.size();
  for(i=0;i<sz;i++) 
    ret.baseCppTypeName_ = ret.baseCppTypeName_ + ret.baseTok_[i] + " ";
  
  ret.wasInstantiated_=true;
  
  return ret;
}


void CppType::copy(CppType const& cppt)
{
  clear_();
  
  cppTypeName_=cppt.cppTypeName_;
  symbolicCppTypeName_=symbolicCppTypeName_;
  baseCppTypeName_=baseCppTypeName_;
  symbolicBaseCppTypeName_=symbolicBaseCppTypeName_;
  
  isConst_=cppt.isConst_;
  isStatic_=cppt.isStatic_;
  isAPointer_=cppt.isAPointer_;
  //  isASTLContainer_=cppt.isASTLContainer_;
  isAComplexType_=cppt.isAComplexType_;
  isAReference_=cppt.isAReference_;
  wasInstantiated_=cppt.wasInstantiated_;
  std::copy(cppt.templateArgs_.begin(),
	    cppt.templateArgs_.end(),
	    back_inserter(templateArgs_));
  std::copy(cppt.tok_.begin(),
	    cppt.tok_.end(),
	    back_inserter(tok_));
  std::copy(cppt.baseTok_.begin(),
	    cppt.baseTok_.end(),
	    back_inserter(baseTok_));
}


void CppType::readFromCppTypeDecl(string const& cpptype)
{
  clear_();
  
  tok_.clear();
  tokenize(cpptype,tok_,",<>*&[] \t", true);
  
  // first pass:
  // analyze the type modifiers
  // and the template arguments
  int i,sz=tok_.size();
  for(i=0;i<sz;) {
    if(tok_[i]=="const") { isConst_=true; i++; continue; }
    if(tok_[i]=="*") { isAPointer_=true; i++; continue; }
    if(tok_[i]=="&") { isAReference_=true; i++; continue; }
    if(tok_[i]=="static") { isStatic_=true; i++; continue; }
    if(tok_[i]=="[") { isAPointer_=true; i++; continue; } // FIXME: should be isArray_
    if(tok_[i]=="<") {
      templateArgs_.push_back(tok_[++i]);
      i++; continue; }
    i++;
  }
  
  // now, try to build the base
  // remove the whitespaces, empty strings and tabs
  // DO NOT remove the ',' (!)
  vector<string>::iterator it;
  for(it=tok_.begin();it!=tok_.end();) {
    
    // remove the blanks
    if(*it==""       || 
       *it==" "      || 
       *it=="\t"     ) { 
      it = tok_.erase(it); 
      continue; 
    }
    
    if(*it=="const"  || 
       *it=="static" || 
       *it=="*"      || 
       *it=="&"      ) {
      it++;
      continue;
    }
    baseTok_.push_back(*it);
    it++;
  }
  
  sz=tok_.size();
  for(i=0;i<sz;i++) 
    cppTypeName_ = cppTypeName_ + tok_[i] + " ";

  sz=baseTok_.size();
  for(i=0;i<sz;i++) 
    baseCppTypeName_ = baseCppTypeName_ + baseTok_[i] + " ";
}



void CppType::print() const
{
  cout << "*** typeName: " << cppTypeName_
       << " (" << baseCppTypeName_  << ")" << endl;
  if(wasInstantiated_)
    cout << "***  [Instantiated from: " << symbolicCppTypeName_ << "]" << endl;
  if(isConst_) cout << "*** const" << endl;
  if(isStatic_) cout << "*** static" << endl;
  if(isAPointer_) cout << "*** pointer type" << endl;
  if(isAReference_) cout << "*** reference" << endl;
  if(isPersistent())
    cout << " ===> PERSISTENT" << endl;
  else
    cout << " ===> NOT PERSISTENT" << endl;
  
  //  int i;
  //  for(i=0;i<tok_.size();i++)
  //    cout << " tok_[" << i << "]=" << tok_[i] << endl;
  
  //  for(i=0;i<baseTok_.size();i++)
  //    cout << " baseTok_[" << i << "]=" << baseTok_[i] << endl;
}




void CppType::tokenize(string const& str,
		       vector<string>& tokens,
		       string const& del, 
		       bool keep_del)
{
  tokens.clear();
  
  string::size_type str_start=str.find_first_not_of(del);
  string::size_type str_end=str.find_first_of(del);
  
  while(str_start != string::npos) {
    if(str_end!=str_start)
      tokens.push_back(str.substr(str_start, str_end-str_start));
    
    if( keep_del ) {
      if(str_end!=string::npos ) {
	tokens.push_back(str.substr(str_end, 1));
	str_start = str_end+1;
	str_end=str.find_first_of(del, str_start);
      }
      else {
	str_start=str_end;
      }
    }
    else {
      str_start=str.find_first_not_of(del, str_end);
      str_end=str.find_first_of(del, str_start);
    }
  }
}



void CppType::clear_()
{
  cppTypeName_="";
  symbolicCppTypeName_="";
  baseCppTypeName_="";
  symbolicBaseCppTypeName_="";
  isConst_=false;
  isStatic_=false;
  isAPointer_=false; 
  // isASTLContainer_=false;
  isAComplexType_=false; 
  isAReference_=false;
  wasInstantiated_=false;
  
  templateArgs_.clear();
  tok_.clear();
  baseTok_.clear();
}


std::string CppTemplateInstance::realName(std::string const& sym) const
{
  int i = find_(sym_,sym);
  if(i<0) return "";
  return real_[i];
}


std::string CppTemplateInstance::symName(std::string const& real) const
{
  int i = find_(real_,real);
  if(i<0) return "";
  return sym_[i];
}


void CppTemplateInstance::addInstance(std::string const& sym, std::string const& real)
{
  if( find_(sym_,sym)>=0 ) {
    std::cout << "CppTemplateInstance::addInstance(sym=" 
	      << sym << ") symbolic type substitution already specified"
	      << std::endl;
    return;
  }
  sym_.push_back(sym);
  real_.push_back(real);
}




void CppTemplateInstance::clear_()
{
  sym_.clear();
  real_.clear();
}


int CppTemplateInstance::find_(std::vector<std::string> const& v, std::string const& str) const
{
  int i,sz=v.size();
  for(i=0;i<sz;i++)
    if(v[i]==str) return i;
  return -1;
}






    //    if(*it=="[") {
    //      it--;
    //      if(*it!="]") {
    //	cout << "YESSS it=" << *it << endl;
    //	it=tok_.erase(it);
    //	//	baseTok_.push_back(*it);
    //	it++;
    //      }
    //      else {
    //	it++;
    //      }
    //      baseTok_.push_back(*it);/      it++;
    //    continue;
    //    }
    //    else {
