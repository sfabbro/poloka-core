// -*- C++ -*-
// 
// dict.cc 
// 
// 
#include <assert.h>

#include <iostream>
#include <iomanip>
#include <algorithm>

#include <config.h>

#include "dict.h"
#include "xmlexceptions.h"
#include "dict_writer_base.h"
#include "dict_reader_base.h"


dict::dict()
  : version_(0), name_(""),
    realName_(""), symName_(""),
    isTemplate_(false),
    isPersistent_(false)
{
}


dict::~dict()
{
  clear_();
}


bool dict::operator==(dict const& d) const
{
  if(size() != d.size()) return false;
  if(version_ != d.version_) return false;
  if(name_ != d.name_) return false;
  if(realName_ != d.realName_) return false;
  if(symName_ != d.symName_) return false;
  if(isTemplate_ != d.isTemplate_) return false;
  if(isPersistent_ != d.isPersistent_) return false;

  unsigned int i;  
  
  for(i=0;i<baseList_.size();i++) 
    if( baseList_[i]!=d.baseList_[i] ) 
      return false;
  
  for(i=0;i<size();i++) 
    if( memberList_[i].name!=d.memberList_[i].name ||
	memberList_[i].type!=d.memberList_[i].type )
      return false;
  
  for(i=0;i<realTypes_.size();i++)
    if( realTypes_[i]!=d.realTypes_[i] )
      return false;
  
  return true;
}


void dict::write(dict_writer_base& dw) const
{
  dw.write(*this);
  //  dw.newDict(name_, version_, isTemplate_, isPersistent_);
  //  dw.writeBaseList(baseList_);
  //  dw.writeMemberList(memberList_);
  //  dw.flush();
}


void dict::read(dict_reader_base const& dr)
{
  dr.read(*this);
  //  dr.newDict(name_, version_, isTemplate_, isPersistent_);
}


void dict::print(int verbosity) const
{
  if(isTemplate_)
    std::cout << "T ";
  else
    std::cout << "  ";
  
  std::cout.flags(ios::left); std::cout.width(10); 
  std::cout << kind_.c_str();
  std::cout << " " << name_ << "[" 
	    << version_ << "] ";
  
  int i,sz=baseList_.size();
  if(sz!=0) {
    std::cout << "{ ";
    for(i=0;i<sz;i++)
      std::cout << baseList_[i] << " ";
    std::cout << "}";
  }
  std::cout << std::endl;
  
  if(verbosity>0) {
    sz=memberList_.size();
    for(i=0;i<sz;i++) 
      memberList_[i].print();
  }
}

//      std::cout << "   + ";
//      std::cout.flags(ios::left); std::cout.width(25); 
//      std::cout << memberList_[i].type.c_str();
//      std::cout.flags(ios::left); std::cout.width(45); 
//      std::cout << memberList_[i].name.c_str()
//		<< std::endl;


void dict::clear_()
{
  version_=0;
  name_="";
  realName_="";
  symName_="";
  isTemplate_=false;
  isPersistent_=false;
  baseList_.clear();
  memberList_.clear();
  symbolicTypes_.clear();
  realTypes_.clear();
}


dict::member_t::member_t()
  : name(""), type(""),
    isPointer(false), isStatic(false),
    isTemplate(false), complexType(false),
    isPersistent(false), isReference(false)
{
}


void dict::member_t::print() const
{
  if(isPersistent)
    std::cout << "   + ";
  else
    std::cout << "   - ";
  std::cout.flags(ios::left); std::cout.width(30); 
  std::cout << type.c_str();
  
  std::cout.flags(ios::left);// std::cout.width(45); 
  std::cout << name.c_str();
  int i,sz=arraySize.size();
  for(i=0;i<sz;i++)
    std::cout << "[" << arraySize[i] << "]";
  
  std::cout << std::endl;
}


void dict::member_t::clear()
{
  name="";
  type="";
  arraySize.clear();
  isPointer=isStatic=isTemplate=complexType=false;
  isPersistent=isReference=false;
}


void dict::member_t::copy(dict::member_t const& d)
{
  clear();
  
  name=d.name;
  type=d.type;
  std::copy(d.arraySize.begin(),d.arraySize.end(),
	    std::back_inserter(arraySize));
  isPointer=d.isPointer;
  isStatic=d.isStatic;
  isTemplate=d.isTemplate;
  complexType=d.complexType;
  isPersistent=d.isPersistent;
  isReference=d.isReference;
}


void dict::member_t::update()
{
  if( isPointer || isReference )
    isPersistent=false;
  else
    isPersistent=true;
  
  if( arraySize.size() == 1 )
    isPersistent=true;
}



void dict::update(templateInstantiation const& ti)
{
  if(!isTemplate_) return;
  
  symName_ = ti.fullSymbolicName(); // to avoid the >> 
  
  // copy the list of symbolic types
  symbolicTypes_.clear();
  int i,sz;
  for(i=0;i<ti.nTemplateArgs();i++)
    symbolicTypes_.push_back(ti.templateSymbolicArg(i));
}


void dict::instantiate(templateInstantiation const& ti)
{
  if(!isTemplate_) return;
  
  symName_ = ti.fullSymbolicName(); 
  realName_ = ti.fullRealName(); // to avoid the >> 
  
  int i,j,sz=size();
  std::string sym;
  
  for(i=0;i<sz;i++) {
    sym = memberList_[i].type;
    memberList_[i].type = ti.instantiateType(sym);
  }
  
  // parse also the baselist
  // which may contain template classes
  sz = baseList_.size();
  for(i=0;i<sz;i++) {
    sym = baseList_[i];
    baseList_[i] = ti.instantiateType(sym);
  }
  
  // copy the list of symbolic and real types
  symbolicTypes_.clear();
  realTypes_.clear();
  for(i=0;i<ti.nTemplateArgs();i++) {
    symbolicTypes_.push_back(ti.templateSymbolicArg(i));
    realTypes_.push_back(ti.templateRealArg(i));
  }
}



std::string dict::templateSymbolicArgList(bool first, bool last) const
{
  int i,sz=symbolicTypes_.size();
  string ret;
  
  if(first && sz)
    ret = ",";
  
  for(i=0;i<sz;i++) 
    ret = ret + symbolicTypes_[i] + ",";
  
  if(!last) {
    std::string::size_type pos = ret.find_last_of(",");
    if(pos!=std::string::npos)
      ret = ret.substr(0,pos);
  }
  return ret;
}


std::string dict::templateSymbolicArgListDecl(bool first, bool last) const
{
  int i,sz=symbolicTypes_.size();
  string ret;
  
  if(first && sz)
    ret = ",";
  
  for(i=0;i<sz;i++) 
    ret = ret + "class " + symbolicTypes_[i] + ",";
  
  if(!last) {
    std::string::size_type pos = ret.find_last_of(",");
    if(pos!=std::string::npos)
      ret = ret.substr(0,pos);
  }
  return ret;
}



templateInstantiation::templateInstantiation()
  : templateClassName_("")
{
}


templateInstantiation::~templateInstantiation()
{
  clear_();
}


void templateInstantiation::clear_()
{
  templateClassName_="";
  symbolicTypes_.clear();
  realTypes_.clear();
}


std::string templateInstantiation::fullSymbolicName() const
{
  std::string ret = templateClassName_ + "<";
  int i, n = nTemplateArgs();
  if(n==0) ret = ret = ">";
  
  for(i=0;i<n-1;i++)
    ret = ret + symbolicTypes_[i] + ",";
  ret = ret + symbolicTypes_[n-1] + ">";
  return ret;
}



std::string templateInstantiation::fullRealName() const
{
  std::string ret = templateClassName_ + "<";
  int i, n = nTemplateArgs();
  if(n==0) ret = ret = ">";
  
  for(i=0;i<n-1;i++)
    ret = ret + realTypes_[i] + ",";
  ret = ret + realTypes_[n-1] + ">";
  return ret;
}


bool templateInstantiation::readFromConfigFile(std::string const& buff)
{
  stringstream sstrm(buff.c_str());
  std::string command, className;
  sstrm >> command >> className;
  if(!sstrm.good()) {
    clear_();
    return false;
  }
  
  if(command!="register")
    return false;
  
  templateClassName_ = className;
  
  while( sstrm.good() ) {
    std::string nm, symtype, realtype;
    std::string::size_type pos;
    sstrm >> nm;
    pos = nm.find_first_of("=");
    symtype = nm.substr(0,pos);
    realtype = nm.substr(pos+1,std::string::npos);
    symbolicTypes_.push_back(symtype);
    realTypes_.push_back(realtype);
  }
  
  return true;
}


void templateInstantiation::print() const
{
  std::cout << "templateInstantiation::print(): " << std::endl;
  std::cout << "  " << fullSymbolicName() 
	    << " --> " << fullRealName()
	    << std::endl;
}


bool templateInstantiation::hasSymbolicType(std::string const& type) const
{
  vector<std::string>::const_iterator it;
  it = std::find(symbolicTypes_.begin(), symbolicTypes_.end(), type);
  return it!=symbolicTypes_.end();
}


std::string templateInstantiation::realType(std::string const& symtype) const
{
  int i,sz=symbolicTypes_.size();
  for(i=0;i<sz;i++)
    if(symbolicTypes_[i]==symtype)
      return realTypes_[i];
  return "";
}


std::string templateInstantiation::instantiateType(std::string const& sym) const
{
  std::string tmp_str=sym;
  std::vector<std::string> vec;
  int i,sz;
  
  splitType_(tmp_str,vec);
  sz=vec.size();
  for(i=0;i<sz;i++) {
    if( hasSymbolicType(vec[i]) ) 
      vec[i] = realType(vec[i]);
  }
  std::string ret = rebuildType_(vec);
  
  // remove the trailing blanks
  std::string::size_type pos = ret.find_last_not_of(" \t");
  if(pos!=std::string::npos) ret = ret.substr(0,pos+1);
  
  return ret;
}


void templateInstantiation::splitType_(std::string const& str, std::vector<std::string>& vec) const
{
  std::string tmp;
  std::string::size_type pos_first=0,pos_last=0,len;
  
  while( pos_first!=std::string::npos ) {
    pos_last=str.find_first_of("<,> \t",pos_first);
    if(pos_last==std::string::npos) {
      len=std::string::npos;
      tmp = str.substr(pos_first,len);
      vec.push_back(tmp);
      return;
    }
    
    len=pos_last-pos_first;
    tmp = str.substr(pos_first,len);
    vec.push_back(tmp);
    tmp = str.substr(pos_last,1);
    vec.push_back(tmp);
    pos_first = pos_last+1;
  }
}


std::string templateInstantiation::rebuildType_(std::vector<std::string> const& vec) const
{
  int i,sz=vec.size();
  std::string ret="";
  for(i=0;i<sz;i++)
    ret = ret + vec[i];
  return ret;
}


