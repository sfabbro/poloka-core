// -*- C++ -*-
// $Id: cpptype.cc,v 1.4 2004/03/08 09:20:29 guy Exp $
// \file cpptype.cc
// 
// 

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "cpptype.h"

using namespace std;


CppType::CppType()
  : isConst_(false),isStatic_(false),
    isAPointer_(false), // isASTLContainer_(false),
    isAComplexType_(false), isAReference_(false),
    isTemplate_(false),
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
  
  //  if(!isTemplate()) {
  //    cout << "CppType::instantiate() ERROR"
  //    	 << " class is not template!" 
  //    	 << endl;
  //    return ret;
  //  }
  
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
  
  
  // OK, now that the substitutions are done,
  // we clean and recompute the type strings
  vector<string> tok;
  string str = buildTypeString(ret.tok_);
  tokenize(str,tok,",<>()*&[] \t", true);
  cleanTokens_(tok);
  ret.tok_.clear();
  std::copy(tok.begin(),tok.end(),
	    std::back_inserter(ret.tok_));
  
  tok.clear();
  str = buildTypeString(ret.baseTok_);
  tokenize(str,tok,",<>()*&[] \t", true);
  cleanTokens_(tok);
  ret.baseTok_.clear();
  std::copy(tok.begin(),tok.end(),
	    std::back_inserter(ret.baseTok_));
  
  
  // now, we recompute the instantiated names
  ret.symbolicCppTypeName_=cppTypeName_;
  ret.symbolicBaseCppTypeName_=baseCppTypeName_;
  ret.cppTypeName_=buildTypeString(ret.tok_);
  ret.baseCppTypeName_ = buildTypeString(ret.baseTok_);
  
  ret.wasInstantiated_=true;
  
  return ret;
}


void CppType::copy(CppType const& cppt)
{
  clear();
  
  cppTypeName_=cppt.cppTypeName_;
  symbolicCppTypeName_=symbolicCppTypeName_;
  baseCppTypeName_=baseCppTypeName_;
  symbolicBaseCppTypeName_=symbolicBaseCppTypeName_;
  
  isConst_=cppt.isConst_;
  isStatic_=cppt.isStatic_;
  isAPointer_=cppt.isAPointer_;
  isAComplexType_=cppt.isAComplexType_;
  isAReference_=cppt.isAReference_;
  isTemplate_=cppt.isTemplate_;
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
  clear();
  
  tok_.clear();
  tokenize(cpptype,tok_,",<>()*&[] \t", true);
  cleanTokens_(tok_);
  
  
  // first pass:
  // analyze the type modifiers
  // and the template arguments
  int i, depth=0, sz=tok_.size();
  for(i=0;i<sz;) {
    if(tok_[i]=="const") { isConst_=true; i++; continue; }
    if(tok_[i]=="*") { isAPointer_=true; i++; continue; }
    if(tok_[i]=="&") { isAReference_=true; i++; continue; }
    if(tok_[i]=="static") { isStatic_=true; i++; continue; }
    if(tok_[i]=="[") { isAPointer_=true; i++; continue; } // FIXME: should be isArray_
    if(tok_[i]=="<") { isTemplate_=true; i++; continue; }
    
    //    if(tok_[i]=="<") {
    //      templateArgs_.push_back(tok_[++i]);
    //      i++; continue; }
    i++;
  }

  // now, the template args.
  // it is a little more tricky
  string tpl_arg="";
  for(i=0;i<sz;i++) {
    if(tok_[i]=="<") {
      if(depth>0)
    	tpl_arg = tpl_arg + tok_[i];
      depth++;
      continue;
    }
    
    if(tok_[i]==">") {
      if(depth==1) {
    	templateArgs_.push_back(tpl_arg);
    	tpl_arg="";
      }
      if(depth>1)
    	tpl_arg = tpl_arg + tok_[i];
      depth--;
      continue;
    }
    
    if(tok_[i]==",") {
      if(depth==1) {
    	templateArgs_.push_back(tpl_arg);
    	tpl_arg="";
      }
      if(depth>1) 
    	tpl_arg = tpl_arg + tok_[i];
      continue;
    }
    
    if(depth>=1)
      tpl_arg = tpl_arg + tok_[i];
  }
  
  // and now, we clean the template args type definition
  int j,ssz=templateArgs_.size();
  for(j=0;j<ssz;j++) {
    std::vector<std::string> tok;
    tokenize(templateArgs_[j],tok,",<>()*&[] \t",true);
    cleanTokens_(tok);
    templateArgs_[j] = buildTypeString(tok);
  }
  
  
  // now, try to build the base type
  // remove the whitespaces, empty strings and tabs
  // DO NOT remove the ',' (!)
  vector<string>::iterator it;
  for(it=tok_.begin();it!=tok_.end();) {
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
  
  cppTypeName_ = buildTypeString(tok_);
  baseCppTypeName_ = buildTypeString(baseTok_);
}



void CppType::print() const
{
  if( isPersistent() )
    std::cout << "   + ";
  else
    std::cout << "   - ";
  std::cout.flags(ios::left); std::cout.width(30); 
  std::cout << cppTypeName().c_str();
  std::cout << std::endl;
}




void CppType::tokenize(string const& str,
		       vector<string>& tokens,
		       string const& del, 
		       bool keep_del)
{
  tokens.clear();
  
  string::size_type str_start=0; // str.find_first_not_of(del); (BUG: if string starts with delimiter...)
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



void CppType::clear()
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
  isTemplate_=false;
  wasInstantiated_=false;
  
  templateArgs_.clear();
  tok_.clear();
  baseTok_.clear();
}


// we remove:
// 1. the blanks
// 2. the <( and )>
void  CppType::cleanTokens_(std::vector<std::string>& tok)
{
  std::vector<std::string>::iterator it,itp;
  
  for(it=tok.begin();it!=tok.end();) {
    if( *it==""  ||
        *it==" " ||
	*it=="(" || // potentiel bug: function pointers are defined with parhenthesis
	*it==")" ||
        *it=="\t" ) {
      it = tok.erase(it);
      continue;
    }
    it++;
  }
  
  for(itp=tok.begin(),it=tok.begin();it!=tok.end();) {
    if(*itp=="<" && *it=="(")
      it=tok.erase(it);
    if(*itp==")" && *it==">")
      it=tok.erase(itp);
    itp=it; it++;
  }
}


std::string CppType::buildTypeString(std::vector<std::string> const& tok)
{
  std::string ret=tok[0];
  
  int i,sz=tok.size();
  
  for(i=1;i<sz;i++) {
    if(tok[i]=="const") {
      ret = ret + " " + tok[i] + " ";
      continue;
    }
    if(tok[i]=="static") {
      ret = ret + " " + tok[i] + " ";
      continue;
    }
    if(tok[i]==">" && tok[i-1]==">") {
      ret = ret + " " + tok[i];
      continue;
    }
    ret = ret + tok[i];
  }
  return ret;
}

std::string CppType::reformatTypeString(std::string const& str)
{
  vector<string> tok;
  tokenize(str,tok,",<>()*&[] \t", true);
  cleanTokens_(tok);
  return buildTypeString(tok);
}


std::string CppType::buildCleanTypeString(std::vector<std::string> const& tok)
{
  std::string ret=tok[0];
  
  int i,sz=tok.size();
  
  for(i=1;i<sz;i++) {
    if(tok[i]==">" ||
       tok[i]=="<" ||
       tok[i]=="(" || 
       tok[i]==")" ||
       tok[i]=="[" || 
       tok[i]=="]" ||
       tok[i]=="*" ||
       tok[i]=="&" ||
       tok[i]=="," ) {
      ret = ret + "_";
      continue;
    }
    ret = ret + tok[i];
  }
  return ret;
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


void CppTemplateInstance::print() const
{
  int i,sz=sym_.size();
  for(i=0;i<sz;i++)
    cout << " " << sym_[i] << " <-> " << real_[i]
	 << endl;
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


//    // remove the blanks
//    if(*it==""       || 
//       *it==" "      || 
//       *it=="\t"     ) { 
//      it = tok_.erase(it); 
//      continue; 
//    }
