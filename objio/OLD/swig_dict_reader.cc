// -*- C++ -*-
// 
// file swig_dict_reader.cc
// 
// 
#include <libxml/parser.h>

#include "swig_dict_reader.h"



swig_dict_reader::swig_dict_reader(std::string const& filename)
  : dict_reader_base(filename), iptr_(0), sz_(0),
    doc_(0), ptr_(0), res_(0), ctx_(0)
{
  doc_ = xmlParseFile(filename_.c_str());
  ptr_ = xmlDocGetRootElement(doc_);
  ctx_ = xmlXPathNewContext(doc_);
  res_ = xmlXPathEvalExpression((xmlChar*)"//class|//template",  ctx_);
  sz_ = res_->nodesetval->nodeNr;
}


swig_dict_reader::~swig_dict_reader()
{
  xmlXPathFreeContext(ctx_);
  xmlFreeDoc(doc_);
}


void swig_dict_reader::read_header_(std::string& name, unsigned int& version,
				    std::string& kind, 
				    bool& isTemplate, bool& isPersistent) const
{
  assert(doc_ && ptr_ && ctx_);
  
  xmlNodePtr node = res_->nodesetval->nodeTab[iptr_];
  ctx_->node = node;
  
  xmlXPathObjectPtr obj;
  
  // name
  obj = 
    xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='name']", ctx_);
  if(!obj || xmlXPathNodeSetIsEmpty(obj->nodesetval) ||
     obj->nodesetval->nodeNr!=1) {
    std::cout << "swig_dict_reader::read_header_ ERROR reading class name"
	      << std::endl;
    return;
  }
  xmlChar* nm = xmlGetProp(obj->nodesetval->nodeTab[0], (xmlChar*)"value");
  name = (const char*)nm;
  xmlFree(nm);
  xmlXPathFreeObject(obj);
  
  // kind
  obj = 
    xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='kind']", ctx_);
  if(!obj || xmlXPathNodeSetIsEmpty(obj->nodesetval)) {
    std::cout << "swig_dict_reader::read_header_ ERROR reading class name"
	      << std::endl;
    return;
  }
  nm = xmlGetProp(obj->nodesetval->nodeTab[0], (xmlChar*)"value");
  kind = (const char*)nm;
  xmlFree(nm);
  xmlXPathFreeObject(obj);
  

  // version
  obj = 
    xmlXPathEvalExpression((xmlChar*)"cdecl/attributelist/attribute[@value='__version__']/..", ctx_);
  if(!obj || xmlXPathNodeSetIsEmpty(obj->nodesetval) ||
     obj->nodesetval->nodeNr!=1) {
    isPersistent=0;
    version=0;
    return;
  }
  xmlNodePtr tmp_node = ctx_->node;
  ctx_->node = obj->nodesetval->nodeTab[0];
  xmlXPathObjectPtr obj_tmp = 
    xmlXPathEvalExpression((xmlChar*)"attribute[@name='value']", ctx_);
  if(!obj_tmp || xmlXPathNodeSetIsEmpty(obj_tmp->nodesetval)) {
    std::cout << "swig_dict::read_header_ ERROR reading class version"
	      << std::endl;
    return;
  }
  nm = xmlGetProp(obj_tmp->nodesetval->nodeTab[0], (xmlChar*)"value");
  version = atoi((const char*)nm);
  xmlFree(nm);
  xmlXPathFreeObject(obj_tmp);
  xmlXPathFreeObject(obj);
  ctx_->node = tmp_node;
  
  if(version>0) isPersistent=true;
  
  std::cout << " version=" << version << std::endl;
}


void swig_dict_reader::read_baselist_(std::vector<std::string>& bl) const
{
  assert(doc_ && ptr_ && ctx_);
  bl.clear();
  
  ctx_->node = res_->nodesetval->nodeTab[iptr_];
  xmlXPathObjectPtr obj;
  
  obj = 
    xmlXPathEvalExpression((xmlChar*)"attributelist/baselist/base", ctx_);
  if(!obj || xmlXPathNodeSetIsEmpty(obj->nodesetval)) {
    return;
  }
  
  int i,sz=obj->nodesetval->nodeNr;
  for(i=0;i<sz;i++) {
    xmlChar* nm = xmlGetProp(obj->nodesetval->nodeTab[i], (xmlChar*)"name");
    bl.push_back((const char*)nm);
    cout << "  bl: " << (const char*)nm << std::endl;
    xmlFree(nm);
  }
    
  
}




void swig_dict_reader::read_memberList_(std::vector<dict::member>& ml) const
{
  assert(doc_ && ptr_ && ctx_);
  ml.clear();  
  
  ctx_->node = res_->nodesetval->nodeTab[iptr_];
  xmlXPathObjectPtr obj, obj2;
  xmlChar* val;
  
  obj = 
    xmlXPathEvalExpression((xmlChar*)"cdecl", ctx_);
  if(!obj) {
    std::cout << "pb reading memberlist."
	      << std::endl;
    return;
  }
  
  
  int i,sz=obj->nodesetval->nodeNr;
  xmlNodePtr p;
  dict::member dm;
  for(i=0;i<sz;i++) {
    dm.name=""; dm.type="";
    
    p = obj->nodesetval->nodeTab[i];
    ctx_->node=p;
    
    obj2 = 
      xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='code']", ctx_);
    if(!obj2 || !xmlXPathNodeSetIsEmpty(obj2->nodesetval)) {
      xmlXPathFreeObject(obj2);
      continue;
    }
    xmlXPathFreeObject(obj2);
    
    obj2 = 
      xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@value='__version__']", ctx_);
    if(!obj2 || !xmlXPathNodeSetIsEmpty(obj2->nodesetval)) {
      xmlXPathFreeObject(obj2);
      continue;
    }
    xmlXPathFreeObject(obj2);
    
    obj2 = 
      xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='name']", ctx_);
    if(!obj2 || xmlXPathNodeSetIsEmpty(obj2->nodesetval)) {
      std::cout << "ERROR parsing member list" << endl;
      return;
    }
    xmlChar* val = xmlGetProp(obj2->nodesetval->nodeTab[0], (xmlChar*)"value");
    dm.name = (const char*)val;
    xmlFree(val);
    xmlXPathFreeObject(obj2);
    
    obj2 = 
      xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='type']", ctx_);
    if(!obj2 || xmlXPathNodeSetIsEmpty(obj2->nodesetval)) {
      std::cout << "ERROR parsing member list" << endl;
      return;
    }
    val = xmlGetProp(obj2->nodesetval->nodeTab[0], (xmlChar*)"value");
    dm.type = (const char*)val;
    xmlFree(val);
    xmlXPathFreeObject(obj2);
    
    if(dm.name!="" && dm.type!="")
      ml.push_back(dm);
    cout << " --> " << dm.type << " " << dm.name << endl;
  }
  
  xmlXPathFreeObject(obj);
  
  cout << " ml.size()=" << ml.size() << endl;
}

//void swig_dict_reader::read(dict& d) const
//{
//  if(iptr_>=sz_) return;
//  std::cout << "read: iptr_=" << iptr_ << std::endl;
//  
//  xmlNodePtr p, node = res_->nodesetval->nodeTab[iptr_];
//  ctx_->node = node;
//  xmlXPathObject res;
//  
//  // class name
//  res = xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='kind']", ctx_);
//  if( xmlXPathNodeSetIsEmpty(res->nodesetval) ) {
//    std::cout << "swig_text_reader_read: ERROR reading class name"
//	      << std::endl;
//    return;
//  }
//  
//  
//  std::cout << " ---> res  : " << res->nodesetval->nodeNr << std::endl;
//  xmlXPathObjectPtr res2 = xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='name']", ctx_);
//  std::cout << " ---> res2 : " << res2->nodesetval->nodeNr << std::endl;
//  xmlXPathObjectPtr res3 = xmlXPathEvalExpression((xmlChar*)"attributelist/baselist/base", ctx_);
//  if(res3 && !xmlXPathNodeSetIsEmpty(res3->nodesetval)) std::cout << " ---> res3 : " << res3->nodesetval->nodeNr << std::endl;
//  xmlXPathObjectPtr res4 = xmlXPathEvalExpression((xmlChar*)"cdecl", ctx_);
//  if(res4 && !xmlXPathNodeSetIsEmpty(res4->nodesetval)) std::cout << " ---> res4 : " << res4->nodesetval->nodeNr << std::endl;
//  
//  iptr_++;
//  std::cout << " returning: iptr_=" << iptr_ << std::endl;
//}

