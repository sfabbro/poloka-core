// -*- C++ -*-
// 
// file xmlschema_dict_writer.cc
// 
// 
#include <config.h>

#include "dict.h"
#include "xmlschema_dict_writer.h"


xmlschema_dict_writer::xmlschema_dict_writer(std::string const& filename)
  : dict_writer_base(filename)
{
  if(filename!="") {
    doc_ = xmlNewDoc((xmlChar*)"1.0");
    rootNode_ = xmlNewNode(0, (xmlChar*)"xsd:schema");
    xsdNs_ = xmlNewNs(rootNode_, (xmlChar*)"http://www.w3.org/2001/XMLSchema", (xmlChar*)"xsd");
    xsiNs_ = xmlNewNs(rootNode_, (xmlChar*)"http://www.w3.org/2001/XMLSchema-instance", (xmlChar*)"xsi");
    xmlDocSetRootElement(doc_,rootNode_);
  }
}


xmlschema_dict_writer::~xmlschema_dict_writer()
{
  if(filename_!="") {
    xmlSaveFormatFileEnc(filename_.c_str(), doc_, "UTF-8", 1);
    xmlFreeDoc(doc_);
  }
}


const char* xsd_comment = 
"\n    XML Schema definition for class %s version %d\n    Written by %s version %s\n      ";

void xmlschema_dict_writer::write(dict const& d)
{
  char buff[1024];
  
  xmlNodePtr cplxTNode = xmlNewChild(rootNode_, xsdNs_,
				     (xmlChar*)"complexType", 0);
  sprintf(buff, "%s", d.name().c_str());
  xmlNewProp(cplxTNode, (xmlChar*)"name", (xmlChar*)buff);
  
  // header and baseList 
  xmlNodePtr annNode = xmlNewChild(cplxTNode, xsdNs_, 
				   (xmlChar*)"annotation", 0);
  sprintf(buff, xsd_comment, d.name().c_str(), d.version(), PACKAGE, VERSION);
  xmlNodePtr docNode = xmlNewChild(annNode, xsdNs_,
				   (xmlChar*)"documentation", (xmlChar*)buff);
  xmlNodePtr appNode = xmlNewChild(annNode, xsdNs_,
				   (xmlChar*)"appinfo", 0);
  sprintf(buff, "%d", d.version());
  xmlNodePtr verNode = xmlNewChild(appNode, 0, (xmlChar*)"version", (xmlChar*)buff);
  sprintf(buff, "%s", d.kind().c_str());
  xmlNodePtr kNode = xmlNewChild(appNode, 0, (xmlChar*)"kind", (xmlChar*)buff);
  xmlNodePtr blNode = xmlNewChild(appNode, 0, (xmlChar*)"baselist", 0);
  int i, sz = d.baseList().size();
  for(i=0;i<sz;i++) {
    sprintf(buff, "%s", d.baseList()[i].c_str());
    xmlNodePtr bNode = xmlNewChild(blNode, 0, (xmlChar*)"base", (xmlChar*)buff);
  }
  
  
  xmlNodePtr cplxCNode = xmlNewChild(cplxTNode, xsdNs_,
				     (xmlChar*)"complexContent", 0);
  xmlNodePtr extNode = xmlNewChild(cplxCNode, xsdNs_,
				   (xmlChar*)"extension", 0);
  xmlNewProp(extNode, (xmlChar*)"base", (xmlChar*)"abstractObjectType");
  xmlNodePtr seqNode = xmlNewChild(extNode, xsdNs_,
				   (xmlChar*)"sequence", 0);
  
  // members 
  sz = d.size();
  for(i=0;i<sz;i++) {
    xmlNodePtr node = xmlNewChild(seqNode, xsdNs_, (xmlChar*)"element", 0);
    xmlNewProp(node, (xmlChar*)"name", (xmlChar*)d.memberName(i).c_str());
    node = xmlNewChild(node, xsdNs_, (xmlChar*)"complexType", 0);
    node = xmlNewChild(node, xsdNs_, (xmlChar*)"complexContent", 0);
    node = xmlNewChild(node, xsdNs_, (xmlChar*)"extension", 0);
    xmlNewProp(node, (xmlChar*)"base", (xmlChar*)d.memberName(i).c_str());
    node = xmlNewChild(node, xsdNs_, (xmlChar*)"attribute", 0);
    xmlNewProp(node, (xmlChar*)"name", (xmlChar*)"name");
    xmlNewProp(node, (xmlChar*)"type", (xmlChar*)"xsd:string");
    xmlNewProp(node, (xmlChar*)"fixed", (xmlChar*)d.memberName(i).c_str());
  }
  
}
