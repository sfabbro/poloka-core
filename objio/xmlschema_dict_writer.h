// -*- C++ -*-
// 
// file xmlschema_dict_writer.h
// 
// 
#ifndef XMLSCHEMA_DICT_WRITER_H
#define XMLSCHEMA_DICT_WRITER_H

#include <string>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include "dict_writer_base.h"


class xmlschema_dict_writer : public dict_writer_base {
public:
  xmlschema_dict_writer(std::string const& filename);
  ~xmlschema_dict_writer();
  
  virtual void  write(dict const&);

private:
  xmlDocPtr doc_;
  xmlNodePtr rootNode_;
  xmlNsPtr xsdNs_;
  xmlNsPtr xsiNs_;
};


#endif

