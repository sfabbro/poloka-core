// -*- C++ -*-
// 
// dict.cc 
// 
// 
#include <assert.h>

#include <iostream>
#include <iomanip>

#include <config.h>

#include "dict.h"
#include "xmlexceptions.h"
#include "dict_writer_base.h"
#include "dict_reader_base.h"


dict::dict()
  : version_(0), name_(""),
    isTemplate_(false), isPersistent_(false)
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
  isTemplate_=false;
  isPersistent_=false;
  baseList_.clear();
  memberList_.clear();
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
  if( isPointer || isReference || arraySize.size() )
    isPersistent=false;
  else
    isPersistent=true;
}

