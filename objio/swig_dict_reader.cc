// -*- C++ -*-
// 
// file swig_dict_reader.cc
// 
// 

#include <iostream>
#include <sstream>
#include <cassert>
#include <libxml/parser.h>

#include "dict.h"
#include "swig_dict_reader.h"

using namespace std;


dict::member_t process_swig_type(std::string const& type);
void tokenize(std::string const& type, 
	      std::vector<std::string>& tokens,
	      std::string const& del=".");
void swig2cpp(std::string& type, dict::member_t& d);



swig_dict_reader::swig_dict_reader(std::string const& filename)
  : dict_reader_base(filename), iptr_(0), sz_(0),
    doc_(0), ptr_(0), res_(0), ctx_(0)
{
  doc_ = xmlParseFile(filename_.c_str());
  ptr_ = xmlDocGetRootElement(doc_);
  ctx_ = xmlXPathNewContext(doc_);
  res_ = xmlXPathEvalExpression((xmlChar*)"/top/include/include/class|/top/include/include/template",  ctx_);
  sz_ = res_->nodesetval->nodeNr;
}



swig_dict_reader::~swig_dict_reader()
{
  //  xmlXPathFreeContext(ctx_);
  //  xmlFreeDoc(doc_);
}


void swig_dict_reader::close()
{
  xmlCleanupParser();
  xmlXPathFreeContext(ctx_); ctx_=0;
  xmlFreeDoc(doc_); doc_=0;
}


void swig_dict_reader::read_header_(std::string& name, unsigned int& version,
				    std::string& kind, 
				    bool& isTemplate, bool& isPersistent) const
{
  assert(doc_ && ptr_ && ctx_);
  
  xmlNodePtr node = res_->nodesetval->nodeTab[iptr_];
  ctx_->node = node;
  
  xmlXPathObjectPtr obj;
  xmlChar* nm;
  
  // template or class/struct ?
  if( xmlStrcmp(res_->nodesetval->nodeTab[iptr_]->name,(xmlChar*)"template")==0 )
    isTemplate=true;
  
  // name
  obj = 
    xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='name']", ctx_);
  //  assert_xpath_selection_(obj,1);
  if(!obj || xmlXPathNodeSetIsEmpty(obj->nodesetval) ||
     obj->nodesetval->nodeNr!=1) {
    std::cout << "swig_dict_reader::read_header_ ERROR reading class name"
  	      << std::endl;
    return;
  }
  nm = xmlGetProp(obj->nodesetval->nodeTab[0], (xmlChar*)"value");
  name = (const char*)nm;
  xmlFree(nm);
  xmlXPathFreeObject(obj);
  
  // kind
  obj = 
    xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='kind']", ctx_);
  //  assert_xpath_selection_(obj,1);
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
}


// read the list of classes we are inheriting from. 
// swig gives us the full type specification, as it appears in the C++ code.
// we have to parse it a little to fill out a dict::member structure.
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
    xmlFree(nm);
  }
}


// read the list of persistent members.
// swig encodes the type its own way.
// again, we have to parse the type in order to fill out 
// a dict::member structure. Remember that if the type is a raw pointer,
// a reference or an array, it won't be marked as persistent (!)
void swig_dict_reader::read_memberList_(std::vector<dict::member_t>& ml) const
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
  std::string swig_type_str;
  
  xmlNodePtr p;
  dict::member_t dm, dm_sw;
  for(i=0;i<sz;i++) {
    dm.name=""; dm.type="";
    
    p = obj->nodesetval->nodeTab[i];
    ctx_->node=p;
    
    obj2 = 
      xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='decl']", ctx_);
    if(!obj2 || xmlXPathNodeSetIsEmpty(obj2->nodesetval)) {
      std::cout << "ERROR parsing member list" << endl;
      return;
    }
    val = xmlGetProp(obj2->nodesetval->nodeTab[0], (xmlChar*)"value");
    std::string decl = (const char*)val;
    xmlFree(val);
    xmlXPathFreeObject(obj2);
    if(decl.find("f(")!=std::string::npos)
      continue;
    swig_type_str = decl;
      
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
    swig_type_str = swig_type_str + (const char*)val;
    xmlFree(val);
    xmlXPathFreeObject(obj2);
    
    dm_sw = process_swig_type(swig_type_str);
    dm_sw.name = dm.name;
    dm_sw.update();
    
    if(dm_sw.name!="" && dm_sw.type!="")
      ml.push_back(dm_sw);
  }
  
  xmlXPathFreeObject(obj);
}



// small utility function to split a string.
// into tokens --like strtok() but safer and easier to use.
// It is used to parse the type specifications.
void tokenize(std::string const& type, 
	      std::vector<std::string>& tokens,
	      std::string const& del)
{
  tokens.clear();
  
  std::string::size_type str_start=type.find_first_not_of(del);
  std::string::size_type str_end=type.find_first_of(del);
  
  while(str_start != std::string::npos) {
    tokens.push_back(type.substr(str_start, str_end-str_start));
    str_start=type.find_first_not_of(del, str_end);
    str_end=type.find_first_of(del, str_start);
  }
}



// convert swig type encoding into C++ type encoding.
void swig2cpp(std::string& type, dict::member_t& d)
{
  // get rid of the ()s
  std::string::size_type pos=0;
  while(pos!=std::string::npos) {
    pos = type.find_first_of("()");
    if(pos!=std::string::npos)
      type.erase(pos,1);
  }
  
  if(type[0]=='q')
    type.erase(0,1);
  else if(type[0]=='a') {
    type.erase(0,1);
    stringstream sstrm(type.c_str());
    unsigned int sz;
    sstrm >> sz;
    d.arraySize.push_back(sz);
    type="";
  }
  else if(type[0]=='p') {
    d.isPointer=true;
    type="*";
  }
  else if(type[0]=='r') {
    d.isReference=true;
    type="&";
  }
  
  d.type = type + " " + d.type;
}


dict::member_t process_swig_type(std::string const& type)
{
  dict::member_t ret;
  
  std::vector<std::string> vv;
  
  tokenize(type, vv);
  int i, sz=vv.size();
  for(i=0;i<sz;i++)
    swig2cpp(vv[i],ret);
  
  // and we remove the trailing blanks
  std::string::size_type pos = ret.type.find_last_not_of(" \t");
  if(pos!=std::string::npos) ret.type = ret.type.substr(0,pos+1);
  
  return ret;
}

////////////////////////////////////////////////////////////////////////////////














//    obj2 = 
//      xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='code']", ctx_);
//    if(!obj2 || !xmlXPathNodeSetIsEmpty(obj2->nodesetval)) {
//      xmlXPathFreeObject(obj2);
//      continue;
//    }
//    xmlXPathFreeObject(obj2);
//    
//    obj2 = 
//      xmlXPathEvalExpression((xmlChar*)"attributelist/parmlist", ctx_);
//    if(!obj2 || !xmlXPathNodeSetIsEmpty(obj2->nodesetval)) {
//      xmlXPathFreeObject(obj2);
//      continue;
//    }
//    xmlXPathFreeObject(obj2);


//void clean_str(std::vector<std::string>& tokens)
//{
//  
//  std::string del = "()";
//  std::vector<std::string>::iterator it;
//  std::string::size_type pos;
//  
//  for(it=tokens.begin();it!=tokens.end();it++) {
//    pos = 0;
//    while(pos!=std::string::npos) {
//      pos = it->find("q(", pos);
//      if(pos!=std::string::npos)
//	it->erase(pos, pos+2);
//    }
    
//    pos = 0;
//    while(pos!=std::string::npos) {
//      pos = it->find("a(", pos);
//      it->erase(pos,pos+1);
//      it->replace(pos,1,"[");
//      pos = it->find(")", pos);
//      it->replace(pos,1);
//    }
    
//    pos = 0;
//    while(pos!=std::string::npos) {
//      pos = it->find_first_of(del);
//      if(pos!=std::string::npos) 
//	it->erase(pos,1);
//    }
//  }
//}






#define SW2CPP(from,to,sz) \
  pos = ret.find(from, 0);\
  if(pos!=std::string::npos) { \
    ret.erase(pos,sz); \
    ret.insert(pos,to); \
    return ret; \
  } \


std::string swig2portablecpp_(std::string const& type)
{
  if(type.size()==1 && 
     type.find("p", 0)!=std::string::npos)
    return "*";
  
  if(type.size()==1 && 
     type.find("r", 0)!=std::string::npos)
    return "&";
  
  //  std::string::size_type pos;
  //  std::string ret=type;
  //  
  //  SW2CPP("unsignedchar","uint1",12);
  //  SW2CPP("unsignedshort","uint2",13);
  //  SW2CPP("unsignedint","uint4",11);
  //  SW2CPP("unsignedlong","uint8",12);
  //  SW2CPP("char","int1",4);
  //  SW2CPP("short","int2",5);
  //  SW2CPP("int","int4",3);
  //  SW2CPP("long","int8",4);
  //  SW2CPP("float","float4",5);
  //  SW2CPP("double","float8",6);
  
  return type;
}


