// -*- C++ -*-
// $Id: cppswigclassreader.h,v 1.2 2004/03/09 10:36:27 nrl Exp $
// \file cppswigclassreader.h
// 
// 
#ifndef CPPSWIGCLASSREADER_H
#define CPPSWIGCLASSREADER_H


#include <libxml/parser.h>
#include <libxml/xpath.h>

#include "cppclassreader.h"
#include "cppclassmember.h"

class CppSwigClassReader : public CppClassReader {
public:
  //  CppSwigClassReader();
  CppSwigClassReader(std::string const& filename);
  ~CppSwigClassReader();
  
  virtual int       size() const { return sz_; }
  virtual CppClass  read() const;
  
  
private:
  mutable unsigned int iptr_;
  unsigned int sz_;
  xmlDocPtr  doc_;
  xmlNodePtr ptr_;
  xmlXPathObjectPtr res_;
  xmlXPathContextPtr ctx_;
  std::vector<bool> classIsOk_;
  
  virtual void next_() const;
  
  virtual void readHeader_(std::string&, unsigned int&, std::string&,
			   bool&, bool&) const;
  virtual void readBaseList_(std::vector<CppClassMember>&) const;
  virtual void readMemberList_(std::vector<CppClassMember>&) const;
  virtual void checkCppClasses_();
};


#endif

