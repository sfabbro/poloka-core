// -*- C++ -*-
// 
// file dict_reader.cc
// 
// 
#include <iostream>
#include <string>
#include <vector>

#include <libxml/parser.h>
#include <libxml/xpath.h>


struct dict {
  dict() : name(""), isTemplate(false) {}
  ~dict() { base.clear(); member.clear(); }
  
  std::string name;
  bool        isTemplate;
  std::vector<std::string> base;
  std::vector<std::string> member;
  void print() const {
    std::cout << "*****************************************" << std::endl
	      << "** name=" << std::setw(29) << " **" << std::endl
	      << "** base.size()=" << base.size() << std::endl
	      << "** member.size()=" << member.size() << std::endl;
  }
};


int main(int argc, char** argv)
{
  xmlDocPtr doc;
  xmlNodePtr cur;
  
  doc = xmlParseFile("test.xml");
  if(doc==0) {
    std::cout << "dict_reader: ERROR reading test.xml"
	      << std::endl;
    return -1;
  }
  
  cur = xmlDocGetRootElement(doc);
  if(cur==0) {
    std::cout << "dict_reader: unable to get the root element"
	      << std::endl;
    return -1;
  }

  xmlXPathContextPtr ctx;
  xmlXPathObjectPtr res;
  ctx = xmlXPathNewContext(doc);
  res = xmlXPathEvalExpression((xmlChar*)"//class|//template", ctx);
  if( xmlXPathNodeSetIsEmpty(res->nodesetval) ) {
    std::cout << "dict_reader: no classes in file"
	      << std::endl;
  }
  
  // are there any classes in there ?
  xmlXPathFreeContext(ctx);
  std::cout << " " << res->nodesetval->nodeNr << " elements found" 
	    << std::endl;
  
  int i,sz=res->nodesetval->nodeNr;
  
  for(i=0;i<sz;i++) {
    read_dict(res->nodesetval->nodeTab[i]);
  
  
  xmlFreeDoc(doc);
  xmlMemoryDump();
}

