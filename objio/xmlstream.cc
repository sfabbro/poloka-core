// -*- C++ -*-
// 
// file xmlstream.cc
// 
// 
#include <assert.h>

#include <list>

#include <config.h>

#include "xmlstream.h"




void  xmlistream::open(std::string const& filename)
{
  if(reader_) close();
  reader_ = xmlNewTextReaderFilename(filename.c_str());
  nextOpeningTag_(); //doc
}


void  xmlistream::close()
{
  if(reader_) {
    xmlFreeTextReader(reader_);
    reader_=0;
  }
}


void xmlistream::start_dict(std::string& name, std::string& kind, 
			    unsigned int& size, unsigned int& version,
			    std::list<std::string>& baselist) const
{
  assert(reader_);
  
  int ret;
  
  
  ret = nextOpeningTag_((xmlChar*)"xsd:annotation");
  ret = nextOpeningTag_((xmlChar*)"xsd:documentation");
  ret = nextOpeningTag_((xmlChar*)"xsd:appinfo");
  ret = nextOpeningTag_((xmlChar*)"version");

  xmlChar const* version_val = xmlTextReaderConstValue(reader_);
  if(!version_val) 
    throw XMLException( BuildExcMsg("unable to read dict version") );
  version = atoi((const char*)version_val);

  ret = nextOpeningTag_((xmlChar*)"kind");
  xmlChar const* kind_val = xmlTextReaderConstValue(reader_);
  if(!kind_val) 
    throw XMLException( BuildExcMsg("unable to read class kind") );
  kind = (const char*)kind_val;
  
  ret = nextOpeningTag_((xmlChar*)"baselist");
  // parsing the list
  baselist.clear();
  while(1) {
    ret = nextElement_();
    if(ret==2) 
      break;
    xmlChar const* base_val = xmlTextReaderConstValue(reader_);
    if(!base_val)
      throw XMLException( BuildExcMsg("unable to read inheritance list") );
    baselist.push_back((const char*)base_val);
    ret = nextElement_();
    if(ret!=2)
      throw XMLException( BuildExcMsg("unable to read inheritance list") );
  }
  
  ret = nextClosingTag_((xmlChar*)"xsd:appinfo");
  ret = nextClosingTag_((xmlChar*)"xsd:annotation");
  ret = nextOpeningTag_((xmlChar*)"xsd:complexContent");
  ret = nextOpeningTag_((xmlChar*)"xsd:extension");
  ret = nextOpeningTag_((xmlChar*)"xsd:sequence");
}


void xmlistream::end_dict() const
{
}


int xmlistream::read_element(std::string& name, std::string& type) const
{
  int ret;
  
  ret = nextElement_();
  if(ret==2)
    return 0;
  
  xmlChar* name_str = xmlTextReaderGetAttribute(reader_, (xmlChar*)"name");
  if(!name_str)
    throw XMLException( BuildExcMsg("error while reading class persistent member def") );
  name = (const char*)name_str;
  xmlFree(name_str);
  
  nextOpeningTag_((xmlChar*)"xsd:complexType");
  nextOpeningTag_((xmlChar*)"xsd:complexContent");
  nextOpeningTag_((xmlChar*)"xsd:extension");  
  xmlChar* type_str = xmlTextReaderGetAttribute(reader_, (xmlChar*)"base");
  if(!type_str)
    throw XMLException( BuildExcMsg("error while reading class persistent member def") );
  type = (const char*)type_str;
  xmlFree(type_str);
  nextOpeningTag_((xmlChar*)"xsd:attribute");  
  
  nextClosingTag_((xmlChar*)"xsd:attribute");
  nextClosingTag_((xmlChar*)"xsd:extension");
  nextClosingTag_((xmlChar*)"xsd:complexContent");
  nextClosingTag_((xmlChar*)"xsd:complexType");
  nextClosingTag_((xmlChar*)"xsd:element");
  
  return 1;
}


////////////////////////////////////////////////////////////////////////////////
// 
// xmlostream stuff...
// 
////////////////////////////////////////////////////////////////////////////////

void xmlostream::open(const std::string& filename, int compression)
{
  if(writer_) close();
  writer_ = xmlNewTextWriterFilename(filename.c_str(), compression);
  xmlTextWriterStartDocument(writer_, "1.0", "UTF-8", "no");
  
  xmlTextWriterWriteFormatString(writer_, "\n\n\n");  
  
  //  xmlTextWriterWriteFormatComment(writer_, "written by %s version %s\n",
  //				  PACKAGE, VERSION);
  
  xmlTextWriterWriteFormatString(writer_, "\n\n\n");
  xmlTextWriterStartElement(writer_, (xmlChar*)"objects");
  xmlTextWriterWriteAttribute(writer_, (xmlChar*)"xmlns:xsi",
			      (xmlChar*)"http://www.w3.org/2001/XMLSchema-instance");
  xmlTextWriterWriteAttribute(writer_, (xmlChar*)"xmlns:snls",
			      (xmlChar*)"http://supernovae.in2p3.fr/SNLS");
  xmlTextWriterWriteFormatString(writer_, "\n");

}


void xmlostream::close()
{
  if(writer_) {
    xmlTextWriterEndDocument(writer_);
    xmlFreeTextWriter(writer_);
    writer_=0;
  }
}



void openTag_(xmlTextWriterPtr w, xmlChar* tagname, unsigned int& depth, bool endline=true)
{
  unsigned int i;
  if(endline) xmlTextWriterWriteFormatString(w,"\n");
  for(i=0;i<depth;i++) xmlTextWriterWriteFormatString(w, " ");
  xmlTextWriterStartElement(w, tagname);
  depth++;
}


void pad_(xmlTextWriterPtr w, unsigned int depth, bool endline=true)
{
  if(endline) xmlTextWriterWriteFormatString(w, "\n");
  unsigned int i;
  for(i=0;i<depth;i++) xmlTextWriterWriteFormatString(w, " ");
}


void closeTag_(xmlTextWriterPtr w, unsigned int& depth, bool endline=true)
{
  if(depth>0) depth--;
  if(endline) xmlTextWriterWriteFormatString(w,"\n");
  unsigned int i;
  for(i=0;i<depth;i++) xmlTextWriterWriteFormatString(w, " ");
  xmlTextWriterEndElement(w);
}



void xmlostream::start_dict(std::string const& name, /* std::string const& kind, */
			    /* unsigned int size, */ unsigned int version
			    /*, std::list<std::string> const& l*/)
{
  assert(writer_);
  
  xmlTextWriterStartElement(writer_, (xmlChar*)"dict");
  xmlTextWriterWriteAttribute(writer_, (xmlChar*)"class", (xmlChar*)name.c_str());
  xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"version", "%d", version);
  //  xmlTextWriterStartElement(writer_, (xmlChar*)"baselist");
  //  std::list<std::string>::const_iterator it;
  //  for(it=l.begin();it!=l.end();it++)
  //    xmlTextWriterWriteFormatElement(writer_, (xmlChar*)"base", "%d", (xmlChar*)(it->c_str()));
  //  xmlTextWriterEndElement(writer_);
}


void xmlostream::write_member(std::string const& name,
			      std::string const& type)
{
  assert(writer_);
  
  xmlTextWriterStartElement(writer_, (xmlChar*)"element");
  xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"name", "%s", (xmlChar*)name.c_str());
  xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"type", "%s", (xmlChar*)type.c_str());
  xmlTextWriterEndElement(writer_);
}


void xmlostream::end_dict()
{
  xmlTextWriterEndElement(writer_);
}




/*
  void  xmlostream::start_dict(std::string const& name, std::string const& kind, 
  unsigned int size, unsigned int version,
  std::list<std::string> const& l)
  {
  assert(writer_);
  
  xmlTextWriterWriteFormatString(writer_, "\n\n");
  
  unsigned int depth=0;
  openTag_(writer_, (xmlChar*)"xsd:annotation", depth);
  xmlTextWriterWriteFormatString(writer_, "\n");
  openTag_(writer_, (xmlChar*)"xsd:documentation", depth);
  pad_(writer_, depth); 
  xmlTextWriterWriteFormatString(writer_, 
				 "XML Schema definition for class %s version %d",
				 name.c_str(), version);
  pad_(writer_, depth);
  xmlTextWriterWriteFormatString(writer_,
				 "Written by %s version %s", PACKAGE, VERSION);
  closeTag_(writer_, depth); // documentation
  
  xmlTextWriterWriteFormatString(writer_, "\n");
  
  openTag_(writer_, (xmlChar*)"xsd:appinfo", depth);
  
  // class version -- crucial
  openTag_(writer_, (xmlChar*)"version", depth);
  pad_(writer_, depth, false);
  xmlTextWriterWriteFormatString(writer_, "%d", version);
  closeTag_(writer_, depth, false);
  
  // whether it is a class or a struct
  openTag_(writer_, (xmlChar*)"kind", depth);
  pad_(writer_, depth, false);
  xmlTextWriterWriteFormatString(writer_, "%s", kind.c_str());
  closeTag_(writer_, depth, false);
  
  // inheritance list
  openTag_(writer_, (xmlChar*)"baselist", depth);
  std::list<std::string>::const_iterator it;
  for(it=l.begin();it!=l.end();it++) {
    openTag_(writer_, (xmlChar*)"base", depth);
    pad_(writer_, depth, false); xmlTextWriterWriteFormatString(writer_, "%s", it->c_str());
    closeTag_(writer_, depth, false);
  }
  closeTag_(writer_, depth); // baselist
  closeTag_(writer_, depth); // appinfo
  
  xmlTextWriterWriteFormatString(writer_, "\n");  
  closeTag_(writer_, depth); // annotation
  
  xmlTextWriterWriteFormatString(writer_, "\n\n\n");
  
  openTag_(writer_, (xmlChar*)"xsd:complexContent", depth);
  openTag_(writer_, (xmlChar*)"xsd:extension", depth);
  xmlTextWriterWriteAttribute(writer_, (xmlChar*)"base", (xmlChar*)"abstractObjectType");
  openTag_(writer_, (xmlChar*)"xsd:sequence", depth);
}


void  xmlostream::end_dict()
{
  assert(writer_);
  
  xmlTextWriterEndElement(writer_);
  xmlTextWriterWriteFormatString(writer_, "\n");
}




void xmlostream::write_element(std::string const& name,
			       std::string const& type)
{
  assert(writer_);
  
  unsigned int depth=2;
  openTag_(writer_, (xmlChar*)"xsd:element", depth);
  xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"name", "%s", name.c_str());
  openTag_(writer_, (xmlChar*)"xsd:complexType", depth);
  openTag_(writer_, (xmlChar*)"xsd:complexContent", depth);
  openTag_(writer_, (xmlChar*)"xsd:extension", depth);
  xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"base", "%s", type.c_str());
  openTag_(writer_, (xmlChar*)"xsd:attribute", depth);
  xmlTextWriterWriteAttribute(writer_, (xmlChar*)"name", (xmlChar*)"name");
  xmlTextWriterWriteAttribute(writer_, (xmlChar*)"type", (xmlChar*)"xsd:string");
  xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"fixed", "%s", (xmlChar*)name.c_str());
  closeTag_(writer_, depth); // attribute
  closeTag_(writer_, depth); // extension
  closeTag_(writer_, depth); // complexContent
  closeTag_(writer_, depth); // complexType
  closeTag_(writer_, depth); // element
  xmlTextWriterWriteFormatString(writer_, "\n");
}
*/













#ifdef GARBAGE
void  xmlistream::start_dict(std::string& name, std::string& kind, 
			     unsigned int& size, unsigned int& version) const
{
  if(!reader_) {
    std::string message = "no file opened";
    // should throw an exc. here
    return;
  }
  
  int ret;
  // should throw exceptions here
  ret = nextOpeningTag_((xmlChar*)"dict");
  
  if(ret<=0) std::cout << "pb: ret=" << ret << std::endl;

  xmlChar* att_val = xmlTextReaderGetAttribute(reader_, (xmlChar*)"className");
  name = (char*)att_val;
  xmlFree(att_val);
  att_val = xmlTextReaderGetAttribute(reader_, (xmlChar*)"kind");
  kind = (char*)att_val;
  xmlFree(att_val);
  att_val = xmlTextReaderGetAttribute(reader_, (xmlChar*)"size");
  size = (unsigned int)atoi((char*)att_val);
  xmlFree(att_val);
  att_val = xmlTextReaderGetAttribute(reader_, (xmlChar*)"version");
  version = (unsigned int)atoi((char*)att_val);
  xmlFree(att_val);
}


void  xmlistream::end_dict() const
{
  if(!reader_) {
    std::string message = "no file opened";
    // should throw an exc. here
    return;
  }
  // next closing tag ?
  nextClosingTag_((xmlChar*)"dict");
}

#endif // GARBAGE
