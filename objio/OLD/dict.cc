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



dict::dict()
  : version_(0), name_(""), kind_(""),
    isTemplate_(false), isPersistent_(false)
{
}


dict::~dict()
{
  clear();
}


//! TODO: rewrite this using the libxml2 native API.
// I do not think that the xmlreader API is very efficient here...
int dict::writeXMLSchema(const std::string& filename) const
{
  xmlTextWriterPtr writer = xmlNewTextWriterFilename(filename.c_str(),0);
  
  xmlTextWriterStartDocument(writer, "1.0", "UTF-8", "yes");
  xmlTextWriterStartElementNS(writer, (xmlChar*)"xsd", (xmlChar*)"schema", 
			      (xmlChar*)"http://www.w3.org/2001/XMLSchema");
  xmlTextWriterWriteAttribute(writer, (xmlChar*)"xmlns:xsi", 
			      (xmlChar*)"http://www.w3.org/2001/XMLSchema-instance");
  
  unsigned int depth=0;
  pad_(writer, depth++);
  xmlTextWriterStartElement(writer, (xmlChar*)"xsd:annotation");
  pad_(writer, depth++);
  xmlTextWriterStartElement(writer, (xmlChar*)"xsd:documentation");
  pad_(writer, depth);
  xmlTextWriterWriteFormatString(writer, "XML Schema definition for class %s version %d", 
				 name().c_str(), version());
  pad_(writer, depth);
  xmlTextWriterWriteFormatString(writer, "Written by %s version %s", PACKAGE, VERSION);
  pad_(writer, depth);
  xmlTextWriterEndElement(writer);
  
  pad_(writer, depth++);
  xmlTextWriterStartElement(writer, (xmlChar*)"xsd:appinfo");
  pad_(writer, depth++);
  xmlTextWriterStartElement(writer, (xmlChar*)"version");
  xmlTextWriterWriteFormatString(writer, "%d", version());
  xmlTextWriterEndElement(writer);
  xmlTextWriterStartElement(writer, (xmlChar*)"kind");
  xmlTextWriterWriteFormatString(writer, "%s", (xmlChar*)kind().c_str());
  xmlTextWriterEndElement(writer);
  
  std::list<std::string>::const_iterator it;
  for(it=inheritanceList_.begin();it!=inheritanceList_.end();it++) {
    pad_(writer, depth);
    xmlTextWriterStartElement(writer, (xmlChar*)"base");
    xmlTextWriterWriteFormatString(writer, "%s", (xmlChar*)(it->c_str()));
    xmlTextWriterEndElement(writer);
  }
  
  pad_(writer, depth--);
  xmlTextWriterEndElement(writer);
  xmlTextWriterEndElement(writer);
  
  
  pad_(writer, depth++); 
  xmlTextWriterStartElement(writer, (xmlChar*)"xsd:complexType");
  xmlTextWriterWriteFormatAttribute(writer, (xmlChar*)"name", "%s", name().c_str());

  pad_(writer, depth++);
  xmlTextWriterStartElement(writer, (xmlChar*)"xsd:complexContent");
  
  pad_(writer, depth++);
  xmlTextWriterStartElement(writer, (xmlChar*)"xsd:extension");
  xmlTextWriterWriteAttribute(writer, (xmlChar*)"base", (xmlChar*)"abstractObjectType");
  
  pad_(writer, depth++);
  xmlTextWriterStartElement(writer, (xmlChar*)"xsd:sequence");
  
  unsigned int i, j;
  for(i=0;i<size();i++) {
    
    pad_(writer,depth);
    pad_(writer,depth);
    
    xmlTextWriterStartElement(writer, (xmlChar*)"xsd:element");   
    xmlTextWriterWriteFormatAttribute(writer, (xmlChar*)"name", "%s", (xmlChar*)(members_[i].name.c_str()));
    pad_(writer, ++depth);
    
    xmlTextWriterStartElement(writer, (xmlChar*)"xsd:complexType");
    pad_(writer, ++depth);
    
    xmlTextWriterStartElement(writer, (xmlChar*)"xsd:complexContent");
    pad_(writer, ++depth);
    
    xmlTextWriterStartElement(writer, (xmlChar*)"xsd:extension");
    xmlTextWriterWriteFormatAttribute(writer, (xmlChar*)"base", "%s", (xmlChar*)(members_[i].typeName.c_str()));
    pad_(writer, ++depth);
    
    xmlTextWriterStartElement(writer, (xmlChar*)"xsd:attribute");
    xmlTextWriterWriteAttribute(writer, (xmlChar*)"name", (xmlChar*)"name");
    xmlTextWriterWriteAttribute(writer, (xmlChar*)"type", (xmlChar*)"xsd:string");
    xmlTextWriterWriteFormatAttribute(writer, (xmlChar*)"fixed", "%s", (xmlChar*)(members_[i].name.c_str()));
    pad_(writer, ++depth);
    
    for(j=0;j<5;j++) {
      xmlTextWriterEndElement(writer);
      pad_(writer, --depth);
    }
  }
  
  xmlTextWriterEndDocument(writer);
  xmlFreeTextWriter(writer);
}


// FIXME: rewrite this with the DOM and XPath -- nrl 02/2004
int dict::readXMLSchema(const std::string& filename)
{
  clear();
  
  xmlTextReaderPtr reader = xmlNewTextReaderFilename(filename.c_str());
  if(reader!=0) 
    throw XMLException( BuildExcMsg("unable to open " + filename) );
  
  int ret;  
  xmlChar* name;

  try {
    ret = nextOpeningTag_(reader, (xmlChar*)"xsd:schema");
    ret = nextOpeningTag_(reader, (xmlChar*)"xsd:annotation");
    ret = nextOpeningTag_(reader, (xmlChar*)"xsd:documentation");
    ret = nextClosingTag_(reader, (xmlChar*)"xsd:documentation");
    ret = nextOpeningTag_(reader, (xmlChar*)"xsd:appinfo");
    ret = nextOpeningTag_(reader, (xmlChar*)"xsd:version");
  } 
  catch(XMLException e) {
    std::cout << "dict::readXMLSchema ERROR reading the dict header" << std::endl;
    throw;
  }
  
  xmlChar* v = xmlTextReaderValue(reader);
  version_ = atoi((const char*)v);
  
  try {
    nextClosingTag_(reader, (xmlChar*)"xsd:version");
    nextClosingTag_(reader, (xmlChar*)"xsd:appinfo");
    nextClosingTag_(reader, (xmlChar*)"xsd:annotation");
    nextOpeningTag_(reader, (xmlChar*)"xsd:complexType");
  }
  catch( XMLException e) {
    std::cout << "dict::readXMLSchema ERROR reading the dict header" << std::endl;
    throw;
  }
  
  name = xmlTextReaderGetAttribute(reader, (xmlChar*)"name");
  if(name==0) {
    xmlFreeTextReader(reader);
    throw XMLException( BuildExcMsg("readXMLSchema ERROR reading name attribute") );
  }
  name_ = (const char*)name; 
  xmlFree(name);
  
  try {
    nextOpeningTag_(reader, (xmlChar*)"xsd:complexContent");
    nextOpeningTag_(reader, (xmlChar*)"xsd:extension");
    nextOpeningTag_(reader, (xmlChar*)"xsd:sequence");
  } 
  catch(XMLException e) {
    xmlFreeTextReader(reader);
    throw XMLException( BuildExcMsg("readXMLSchema ERROR reading the class body") );
  }
  
  unsigned int i, j;
  while(1) {
    persistentMember pm;
    
    ret = nextOpeningTag_(reader);
    if(ret==0 || ret == 2)
      break;
    
    name = xmlTextReaderGetAttribute(reader, (xmlChar*)"name");
    if(name==0) {
      std::cout << "readXMLSchema ERROR reading name attribute" << std::endl;
      xmlFreeTextReader(reader);
      throw XMLException( BuildExcMsg("unable to read name attribute") );
    }
    pm.name = (const char*)name; 
    xmlFree(name);
    
    try {
      nextOpeningTag_(reader, (xmlChar*)"xsd:complexType");
      nextOpeningTag_(reader, (xmlChar*)"xsd:complexContent");
      nextOpeningTag_(reader, (xmlChar*)"xsd:extension");
    }
    catch( XMLException e ) {
      std::cout << "readXMLSchema ERROR reading element=" 
		<< (const char*)name << std::endl;
      throw;
    }
    
    name = xmlTextReaderGetAttribute(reader, (xmlChar*)"base");
    if(name==0) {
      xmlFreeTextReader(reader);
      throw XMLException( BuildExcMsg("readXMLSchema ERROR reading base attribute") );
    }
    
    pm.typeName = (const char*)name; 
    xmlFree(name);
    addPersistentMember(pm);
    
    try {
      nextOpeningTag_(reader, (xmlChar*)"xsd:attribute");
      nextClosingTag_(reader, (xmlChar*)"xsd:attribute");
      nextClosingTag_(reader, (xmlChar*)"xsd:extension");
      nextClosingTag_(reader, (xmlChar*)"xsd:complexContent");
      nextClosingTag_(reader, (xmlChar*)"xsd:complexType");    
      nextClosingTag_(reader, (xmlChar*)"xsd:element");
    }
    catch( XMLException e ) {
      throw;
    }
  }
  
  xmlFreeTextReader(reader);
}




int dict::readFromSwigXMLStream(xmlTextReaderPtr reader)
{
  clear();
  
  assert(reader);
  
  int ret;
  xmlChar* name;
  
  ret = nextElement_(reader, (xmlChar*)"class");
  if(ret<0) {
    std::cout << "dict::readFromSwigXMLStream ERROR getting to next <class> tag."
	      << std::endl;
    return ret;
  }
  
  if(ret==0)
    return ret;
  
  ret = readClassAttributeList_(reader);

  if(ret<=0) {
    std::cout << "[readFromSwigXMLStream] ERROR parsing the class attribute list."
	      << std::endl;
  }
  
  ret = readClassBody_(reader);
  
  return 1;
}


int dict::nextOpeningTag_(xmlTextReaderPtr reader, xmlChar const* tagname)
{
  assert(reader);
  
  int ret, type;
  while(1) {
    ret = xmlTextReaderRead(reader);
    if(ret==0) 
      return ret;
    if(ret<0)
      throw XMLException( BuildExcMsg("error reading the next opening tag") );
    type = xmlTextReaderNodeType(reader);
    if(type==XML_READER_TYPE_ELEMENT) 
      break;
  }
  if(!tagname) 
    return 1;
  xmlChar const* cur_tagname = xmlTextReaderConstName(reader);
  if( xmlStrcmp(tagname, cur_tagname)==0 ) 
    return 1;
  throw XMLException( BuildExcMsg("format error: unexpected tag") );
}


int dict::nextClosingTag_(xmlTextReaderPtr reader, xmlChar const* tagname)
{
  assert(reader);
  
  int ret, type;
  while(1) {
    xmlTextReaderRead(reader);
    if(ret==0) 
      return ret;
    if(ret<0)
      throw XMLException( BuildExcMsg("error reading the next closing tag") );
    type=xmlTextReaderNodeType(reader);
    if(type==XML_READER_TYPE_END_ELEMENT) 
      break;
  }
  if(!tagname) 
    return 1;
  xmlChar const* cur_tagname = xmlTextReaderConstName(reader);
  if( xmlStrcmp(tagname, cur_tagname) == 0 )
    return 1;
  throw XMLException( BuildExcMsg("format error: unexpected tag") );
}


int dict::nextElement_(xmlTextReaderPtr reader, xmlChar const* name)
{
  assert(reader);
  int ret, type, final_ret;
  
  while(1) {
    ret = xmlTextReaderRead(reader);
    if(ret==0) 
      return ret;
    if(ret<0) 
      throw XMLException("error reading the next element");
    type = xmlTextReaderNodeType(reader);
    if(type == XML_READER_TYPE_ELEMENT ) 
      final_ret=1;
    else if(type == XML_READER_TYPE_END_ELEMENT) 
      final_ret=2;
    else 
      continue;
    if(!name) 
      return final_ret;
    xmlChar const* curName = xmlTextReaderConstName(reader);
    if( xmlStrcmp(curName, name) == 0 ) 
      return final_ret; 
  }
}

int dict::readClassAttributeList_(xmlTextReaderPtr reader)
{
  assert(reader);
  
  int ret;
  std::string name_str, value_str;
  
  ret = nextOpeningTag_(reader, (xmlChar*)"attributelist");
  if(ret<=0) {
    std::cout << "[readClassAttributeList_] ERROR tag <attributelist> not found"
	      << std::endl;
    return ret;
  }
  
  while(1) {
    ret = nextElement_(reader);
    
    if(ret<=0) {
      std::cout << "[readClassAttributeList_] ERROR parsing the class attribute list"
		<< std::endl;
      return -1;
    }

    const xmlChar* name = xmlTextReaderConstName(reader);
    
    if(ret==2) {
      if( xmlStrcmp(name,(xmlChar*)"attributelist") ==0 )
	break;
    }
    
    if( xmlStrcmp(name, (xmlChar*)"attribute") == 0 ) {
      decodeAttribute_(reader, name_str, value_str);
      if(name_str=="name") name_ = value_str;
      if(name_str=="kind") kind_ = value_str;
    }
    else if( xmlStrcmp(name, (xmlChar*)"baselist") == 0 )
      readInheritanceList_(reader);
  }
  
  return 1;
}



int dict::readInheritanceList_(xmlTextReaderPtr reader)
{
  assert(reader);
  
  int ret;
  
  while(1) {
    ret = nextElement_(reader);
    if(ret<=0) {
      std::cout << "[readInheritanceList_] ERROR while parsing class Inheritance list"
		<< std::endl;
      return -1;
    }

    xmlChar const* tagname = xmlTextReaderConstName(reader);
    
    if(ret==2) {
      if( xmlStrcmp(tagname,(xmlChar*)"baselist")==0 )
	break;
    }
    
    xmlChar* attname = xmlTextReaderGetAttribute(reader, (xmlChar*)"name");
    inheritanceList_.push_back((char*)attname);
    xmlFree(attname);
  }
  
  return 1;
}



int dict::readClassBody_(xmlTextReaderPtr reader)
{
  assert(reader);
  
  int ret;
  
  while(1) {
    ret = nextElement_(reader);
    
    if(ret<=0) {
      std::cout << "dict::readClassBody_ ERROR end of file reached."
		<< std::endl;
      return ret;
    }
    
    xmlChar const* tagname=xmlTextReaderConstName(reader);
    
    if(ret==2) {
      if( xmlStrcmp(tagname, (xmlChar*)"class") == 0 )
	break;
      else 
	continue;
    }
    
    if(ret!=1) {
      std::cout << "dict::readClassBody_ ERROR unknown return code ret="
		<< ret 
		<< std::endl;
      return -1;
    }
    
    if( xmlStrcmp(tagname, (xmlChar*)"cdecl")!=0 )
      continue;
    
    readCDeclAttributeList_(reader);
  }
  
  return 1;
}

int dict::readCDeclAttributeList_(xmlTextReaderPtr reader)
{
  assert(reader);
  
  int ret;
  std::string name_str, value_str;
  std::string el_name="", el_type="", el_value="";
  std::string final_type_name="";
  bool skip = false;
  
  ret = nextOpeningTag_(reader, (xmlChar*)"attributelist");
  if(ret<=0) {
    std::cout << "[readCDeclAttributeList_] ERROR tag <attributelist> not found"
	      << std::endl;
    return ret;
  }
  
  
  // FIXME: throw away thse depth tests and 
  // implement a real nextSibling...
  //  int cur_depth, depth=xmlTextReaderDepth(reader);
  
  while(1) {
    ret = nextElement_(reader);
    if(ret<=0) {
      std::cout << "[readCDeclAttributeList_] ERROR parsing CDecl attribute list"
		<< std::endl; 
      return -1;
    }
    
    const xmlChar* tagname = xmlTextReaderConstName(reader);    
    
    if(ret==2)
      if( xmlStrcmp(tagname, (xmlChar*)"attributelist") == 0 )
	break;
    
    if( xmlStrcmp(tagname, (xmlChar*)"attribute") == 0 ) {
      decodeAttribute_(reader, name_str, value_str);
      if(name_str=="name") el_name = value_str;
      if(name_str=="type") el_type = value_str;
      if(name_str=="value") el_value = value_str;
      if(name_str=="code") skip = true;
    }
  }

  if(el_name=="" || el_type=="" || skip) return 1;
  if(el_name=="__version__" && el_type=="q(const).unsigned short") {
    isPersistent_=true;
    version_=(unsigned int)atoi(el_value.c_str());
  }
  else {
    final_type_name = tr_type_name_(el_type);
    if(final_type_name!="UNKNOWN") {
      persistentMember pm;
      pm.name = el_name;
      pm.typeName = final_type_name;
      members_.push_back(pm);
    }
  }
  
  return 1;
}


void dict::decodeAttribute_(xmlTextReaderPtr reader, std::string& name, std::string& value)
{
  xmlChar* name_tmp = xmlTextReaderGetAttribute(reader, (xmlChar*)"name");
  xmlChar* value_tmp = xmlTextReaderGetAttribute(reader, (xmlChar*)"value");
  
  name = (char*)name_tmp;
  value = (char*)value_tmp;
  
  xmlFree(name_tmp);
  xmlFree(value_tmp);
}


bool dict::operator=(dict const& d) const 
{
  if( (name_!=d.name_) ) return false;
  unsigned int i;
  for(i=0;i<members_.size();i++)
    if( !(members_[i]=d.members_[i]) )
      return false;
  return true;
}


void dict::print(int verbosity) const
{
  std::cout << std::resetiosflags(std::ios::left);
  if(isTemplate_)
    std::cout << "T ";
  else
    std::cout << "  ";
  
  std::cout << std::setw(2) << version_ << " ";
  
  std::cout << std::setw(7) << kind_ << " "
	    << name_;
  if(inheritanceList_.size()>0) {
    std::cout << " [" ;
    std::list<std::string>::const_iterator it;
    for(it=inheritanceList_.begin();it!=inheritanceList_.end();it++)
      std::cout << " " << *it ;
    std::cout << " ]" << std::endl;
  }
  else {
    std::cout << std::endl;
  }
  
  if(verbosity>0) {
    unsigned int i;
    for(i=0;i<members_.size();i++) {
      std::cout << "             ";
      members_[i].print();
    }
  }
}



void dict::clear()
{
  version_=0;
  name_="";
  kind_="";
  isTemplate_=false;
  isPersistent_=false;
  inheritanceList_.clear();
  members_.clear();
}


void dict::test(xmlTextReaderPtr reader)
{
  std::cout << " test!" << std::endl;

  int ret = nextElement_(reader, (xmlChar*)"class");
  std::cout << "OK///" << std::endl;
  if(ret<=0) std::cout << " no found..." << std::endl;
  
  while(ret>0) {
    ret = nextOpeningTag_(reader);
    xmlChar const* n1 = xmlTextReaderConstName(reader);
    ret = nextClosingTag_(reader);
    xmlChar const* n2 = xmlTextReaderConstName(reader);
    std::cout << " n1=" << (char*)n1 << ", n2=" << (char*)n2 << std::endl;
  }
}


void dict::printNodeInfo_(xmlTextReaderPtr reader)
{
  int type = xmlTextReaderNodeType(reader);
  xmlChar const* tagname=xmlTextReaderConstName(reader);
  
  switch(type) {
  case XML_READER_TYPE_ELEMENT:
    std::cout << "dict::printNodeInfo() <" << (char*)tagname 
	      << ">" << std::endl;
    break;
  case XML_READER_TYPE_END_ELEMENT:
    std::cout << "dict::printNodeInfo() </" << (char*)tagname 
	      << ">" << std::endl;
    break;
  default:
    std::cout << " UNKNOWN, name=" << (char*)tagname
	      << std::endl;
  }
}




void dict::pad_(xmlTextWriterPtr writer, unsigned int depth) const
{
  unsigned int i;
  xmlTextWriterWriteFormatString(writer, "\n");
  for(i=0;i<depth;i++)
    xmlTextWriterWriteFormatString(writer, " ");
}




std::string dict::tr_type_name_(std::string const& stn)
{
  if(stn=="int1") return (std::string)"i1";
  if(stn=="uint1") return (std::string)"u1";
  if(stn=="int2") return (std::string)"i2";
  if(stn=="uint2") return (std::string)"u2";
  if(stn=="int4") return (std::string)"i4";
  if(stn=="uint4") return (std::string)"u4";
  if(stn=="int8") return (std::string)"i8";
  if(stn=="uint8") return (std::string)"u8";
  if(stn=="float4") return (std::string)"f4";
  if(stn=="float8") return (std::string)"f8";
  if(stn=="string") return (std::string)"string";
  std::cout << "dict::tr_type_name: UNKNOWN type:" <<  stn
	    << " skipped." << std::endl;
  return (std::string)"UNKNOWN";
}




