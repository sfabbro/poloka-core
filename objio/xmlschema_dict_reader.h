// -*- C++ -*-
//
// file xmlschema_dict_reader.h
//
//
#ifndef XMLSCHEMA_DICT_READER_H
#define XMLSCHEMA_DICT_READER_H

#include <libxml/parser.h>
#include <libxml/xpath.h>

#include "dict_reader_base.h"

class xmlschema_dict_reader : public dict_reader_base {
 public:
  xmlschema_dict_reader(std::string const& filename);
  ~xmlschema_dict_reader();

  virtual unsigned int size() const { return sz_; }

 private:
  mutable unsigned int iptr_;
  unsigned int sz_;
  xmlDocPtr  doc_;
  xmlNodePtr ptr_;
  xmlXPathObjectPtr res_;
  xmlXPathContextPtr ctx_;


  virtual void next_() const { if(iptr_>=sz_) return; iptr_++; }
  virtual void read_header_(std::string&, unsigned int&, std::string&,
                            bool&, bool&) const;
  virtual void read_baselist_(std::vector<std::string>&) const;
  virtual void read_memberList_(std::vector<dict::member>&) const;
};


#endif

