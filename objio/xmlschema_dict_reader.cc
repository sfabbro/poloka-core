// -*- C++ -*-
// 
// file xmlschema_dict_reader.cc
// 
// 
#include "xmlschema_dict_reader.h"


xmlschema_dict_reader::xmlschema_dict_reader(std::string const& filename)
  : dict_reader_base(filename), 
    doc_(0), ptr_(0), res_(0), ctx_(0)
{
  doc_ = xmlParserFile(filename_.c_str());
  ptr_ = xmlDocGetRootElement(doc_);
  ctx_ = xmlXPathNewContext(doc_);
  res_ = xmlXPathEvalExpression((xmlChar*)"schema/complexContent", ctx_);
  sz_ = res_->nodesetval->nodeNr;
}


xmlschema_dict_reader::~xmlschema_dict_reader()
{
  xmlXPathFreeContext(ctx_);
  xmlFreeDoc(doc_);
}



void xmlschema_dict_reader::reader_header_(std::string& name, unsigned int& version,
					   std::string& kind,
					   bool& isTemplate, bool& isPersistent)
{
  assert(doc_ && ptr_ && ctx_);
  xmlChar* val;  
  
  xmlNodePtr node = res_->nodesetval->nodeTab[iptr_];
  ctx_->node = node;
  
  xmlGetProp(node, (xmlChar*)"name");
  
  xmlXPathObjectPtr obj;
  
  
  obj = 
    xmlXPathEvalExpression((xmlChar*)"documentation/appinfo/version");
  val = xmlNodeListGetString(doc_, obj->nodesetval->nodeTab[0], 1);
  version = atoi((const char*)val);
  xmlFree(val);
  xmlXPathFreeObject(obj);
  
  obj = 
    xmlXPathEvalExpression((xmlChar*)"kind");
  val = xmlNodeListGetString(doc_, obj->nodesetval->nodeTab[0], 1);
  kind = (const char*)val;
  xmlFree(val);
  xmlXPathFreeObject(obj);
  
}

