// -*- C++ -*-
// $Id: cppswigclassreader.cc,v 1.2 2004/03/08 13:20:01 guy Exp $
// \file cppswigclassreader.cc
// 
// 
#include <iostream>
#include <string>
#include <sstream>

#include "cpptype.h"
#include "cppswigclassreader.h"


using namespace std;

string swig2cpp(string& type);


CppSwigClassReader::CppSwigClassReader(string const& filename)
  : CppClassReader(filename), iptr_(0), sz_(0),
    doc_(0), ptr_(0), res_(0), ctx_(0)
{
  doc_ = xmlParseFile(filename_.c_str());
  ptr_ = xmlDocGetRootElement(doc_);
  ctx_ = xmlXPathNewContext(doc_);
  res_ = xmlXPathEvalExpression((xmlChar*)"/top/include/include/class|/top/include/include/template",  ctx_);
  sz_ = res_->nodesetval->nodeNr;
}


CppSwigClassReader::~CppSwigClassReader()
{
  xmlCleanupParser();
  xmlXPathFreeContext(ctx_); ctx_=0;
  xmlFreeDoc(doc_); doc_=0;
}


CppClass CppSwigClassReader::read() const
{
  int i,sz;
  string name, kind;
  unsigned int version;
  bool isTemplate, isPersistent;
  vector<CppClassMember> base;
  vector<CppClassMember> members;
  
  readHeader_(name, version, kind, isTemplate, isPersistent);
  readBaseList_(base);
  readMemberList_(members);
  
  CppClass ret(kind, name, version, isTemplate);
  sz=base.size();
  for(i=0;i<sz;i++) ret.addBaseClass(base[i]);
  sz=members.size();
  for(i=0;i<sz;i++) ret.addMember(members[i]);
  
  next_();
  
  return ret;
}


void CppSwigClassReader::readHeader_(string& name, unsigned int& version,
				     string& kind,
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
  else
    isTemplate=false;
  
  // name
  obj = 
    xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='name']", ctx_);
  //  assert_xpath_selection_(obj,1);
  if(!obj || xmlXPathNodeSetIsEmpty(obj->nodesetval) ||
     obj->nodesetval->nodeNr!=1) {
    cout << "swig_dict_reader::read_header_ ERROR reading class name"
	 << endl;
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
    cout << "swig_dict_reader::read_header_ ERROR reading class name"
	 << endl;
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
    cout << "swig_dict::read_header_ ERROR reading class version"
	 << endl;
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
void CppSwigClassReader::readBaseList_(vector<CppClassMember>& bl) const
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
    CppClassMember cm((const char*)nm,(const char*)nm);
    bl.push_back(cm);
    xmlFree(nm);
  }
}




// read the list of persistent members.
// swig encodes the type its own way.
// again, we have to parse the type in order to fill out 
// a dict::member structure. Remember that if the type is a raw pointer,
// a reference or an array, it won't be marked as persistent (!)
void CppSwigClassReader::readMemberList_(vector<CppClassMember>& ml) const
{
  assert(doc_ && ptr_ && ctx_);
  ml.clear();  
  
  ctx_->node = res_->nodesetval->nodeTab[iptr_];
  xmlXPathObjectPtr obj, obj2;
  xmlChar* val;
  
  obj = 
    xmlXPathEvalExpression((xmlChar*)"cdecl", ctx_);
  if(!obj) {
    cout << "pb reading memberlist."
	 << endl;
    return;
  }
  
  
  int i,sz=obj->nodesetval->nodeNr;
  string swig_type_str;
  
  xmlNodePtr p;
  string member_name;
  string member_type;
  //  dict::member_t dm, dm_sw;
  
  for(i=0;i<sz;i++) {
    member_name="";
    //    dm.name=""; dm.type="";
    
    p = obj->nodesetval->nodeTab[i];
    ctx_->node=p;
    
    // first, get to the next <cdecl> element
    // and try to determine whether it is a function or a data member.
    obj2 = 
      xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='decl']", ctx_);
    if(!obj2 || xmlXPathNodeSetIsEmpty(obj2->nodesetval)) {
      cout << "ERROR parsing member list" << endl;
      return;
    }
    val = xmlGetProp(obj2->nodesetval->nodeTab[0], (xmlChar*)"value");
    string decl = (const char*)val;
    xmlFree(val);
    xmlXPathFreeObject(obj2);
    if(decl.find("f(")!=string::npos)
      continue;
    swig_type_str = decl;
    
    // if the element name is __version__,
    // it is not a data member, but a field indicating 
    // the class version. It was read already in the readHeader_() method
    // So, we skip it.
    obj2 = 
      xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@value='__version__']", ctx_);
    if(!obj2 || !xmlXPathNodeSetIsEmpty(obj2->nodesetval)) {
      xmlXPathFreeObject(obj2);
      continue;
    }
    xmlXPathFreeObject(obj2);
    
    // now, try to find the member name
    obj2 = 
      xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='name']", ctx_);
    if(!obj2 || xmlXPathNodeSetIsEmpty(obj2->nodesetval)) {
      cout << "ERROR parsing member list" << endl;
      return;
    }
    xmlChar* val = xmlGetProp(obj2->nodesetval->nodeTab[0], (xmlChar*)"value");
    member_name=(const char*)val;
    xmlFree(val);
    xmlXPathFreeObject(obj2);
    
    // and the member type
    obj2 = 
      xmlXPathEvalExpression((xmlChar*)"attributelist/attribute[@name='type']", ctx_);
    if(!obj2 || xmlXPathNodeSetIsEmpty(obj2->nodesetval)) {
      cout << "ERROR parsing member list" << endl;
      return;
    }
    val = xmlGetProp(obj2->nodesetval->nodeTab[0], (xmlChar*)"value");
    swig_type_str = swig_type_str + (const char*)val;
    xmlFree(val);
    xmlXPathFreeObject(obj2);
    
    member_type = swig2cpp(swig_type_str);
    
    // if we have a valid name and a valid type
    // we save the element:
    if(member_type!="" && member_name!="") {
      CppClassMember m(member_type, member_name);
      //      cout << " --> " << member_type << endl;
      //      cout << " ++> "; m.print();
      ml.push_back(m);
    }
    //    if(dm_sw.name!="" && dm_sw.type!="")
    //      ml.push_back(dm_sw);
  }
  
  xmlXPathFreeObject(obj);
}




// convert swig type encoding into C++ type encoding.
string swig2cpp(string& type)
{
  vector<string> tok, array_size;
  
  CppType::tokenize(type,tok,". \t",false);
  string str, ret;
  
  int i,sz=tok.size();
  
  
  for(i=0;i<sz;i++) {
    str=tok[i];
      
    // get rid of the ()s
    string::size_type pos=0;
    while(pos!=string::npos) {
      pos = str.find_first_of("()");
      if(pos!=string::npos)
	str.erase(pos,1);
    }
    
    if(str[0]=='q')
      str.erase(0,1);
    else if(str[0]=='a') {
      str.erase(0,1);
      array_size.push_back(str);
      continue;
    }
    else if(str[0]=='p') {
      str="*";
    }
    else if(str[0]=='r') {
      str="&";
    }
    
    ret = ret + " " + str;
  }
  
  sz = array_size.size();
  for(i=0;i<sz;i++)
    ret = ret + "[" + array_size[i] + "]";
  
  return ret;
}



