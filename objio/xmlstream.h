// -*- C++ -*-
// 
// $Id: xmlstream.h,v 1.3 2004/02/27 16:34:36 nrl Exp $
// 
// 
#ifndef XMLSTREAM_H
#define XMLSTREAM_H

#include <assert.h>

#include <iostream>
#include <list>
#include <string>

#include <libxml/encoding.h>
#include <libxml/xmlreader.h>
#include <libxml/xmlwriter.h>

#include "xmlexceptions.h"
#include "toadtypes.h"




class xmlistream {
public:
  xmlistream() : reader_(0) {}
  xmlistream(std::string const& filename) 
    : reader_(0) { open(filename); }
  ~xmlistream() { close(); }
  
  void            open(std::string const&);
  void            close();
  
  inline void     start_object(std::string& name, unsigned int& version) const;
  inline void     end_object() const;
  
  inline void     start_raw_pointer() const;
  inline void     end_raw_pointer() const;

  inline void     start_reference() const;
  inline void     end_reference() const;
  
  inline void     start_collection(unsigned int&) const;
  inline void     end_collection() const;
  
  inline void     read(int1& v) const;
  inline void     read(uint1& v) const;
  inline void     read(int2& v) const;
  inline void     read(uint2& v) const;
  inline void     read(int4& v) const;
  inline void     read(uint4& v) const;
  inline void     read(int8& v) const;
  inline void     read(uint8& v) const;
  inline void     read(float4& v) const;
  inline void     read(float8& v) const;
  inline void     read(std::string&) const;
  
  inline void     skip() const {
    nextOpeningTag_();
    nextClosingTag_();
  }
  
  void            start_dict(std::string& name, std::string& kind, 
			     unsigned int& size, unsigned int& version,
			     std::list<std::string>& baselist) const;
  
  void            end_dict() const;
  
  int             read_element(std::string& name, std::string& type) const;
  
private:
  xmlTextReaderPtr reader_;
  
  inline int nextOpeningTag_(xmlChar const* tagname=0) const throw(XMLException);
  inline int nextClosingTag_(xmlChar const* tagname=0) const throw(XMLException);
  inline int nextTextElement_() const throw(XMLException);
  inline int nextElement_(xmlChar const* tagname=0) const throw(XMLException);
  
  void readElementName_(std::string& name) const {
    xmlChar const* nm = xmlTextReaderConstName(reader_);
    name = (char*)nm;
  }
  
  void readElementVersion_(unsigned int& version) const {
    version=0;
    xmlChar* ver_str = xmlTextReaderGetAttribute(reader_, (xmlChar*)"version");
    if(ver_str==0) return; // should throw an exc. here
    version = (unsigned int)atoi((char*)ver_str);
  }
};





class xmlostream {
public:
  xmlostream() : writer_(0) {}
  xmlostream(std::string const& filename, int compression=0) 
    : writer_(0) { open(filename, compression); }
  ~xmlostream() { close(); }
  
  void            open(std::string const& filename, int compression=0);
  void            close();
  
  inline void     start_object(std::string const& name, unsigned int version);
  inline void     end_object();
  
  inline void     start_raw_pointer(const char* name=0);
  inline void     end_raw_pointer();
  
  inline void     start_reference(const char* name=0);
  inline void     end_reference();
  
  inline void     start_collection(unsigned int size, const char* name=0);
  inline void     end_collection();
  
  inline void     write(int1 v, const char* name=0);
  inline void     write(uint1 v, const char* name=0);
  inline void     write(int2 v, const char* name=0);
  inline void     write(uint2 v, const char* name=0);
  inline void     write(int4 v, const char* name=0);
  inline void     write(uint4 v, const char* name=0);
  inline void     write(int8 v, const char* name=0);
  inline void     write(uint8 v, const char* name=0);
  inline void     write(float4 v, const char* name=0);
  inline void     write(float8 v, const char* name=0);
  inline void     write(const std::string& v, const char* name=0);
  
  // meta data
  void            start_dict(std::string const& name, /* std::string const& kind,*/
			     /* unsigned int size, */ unsigned int version
			     /*, std::list<std::string> const& baselist */);
  void            end_dict();
  
  void            write_member(std::string const& name,
			       std::string const& type);
  
private:
  xmlTextWriterPtr writer_;
  
  void   writeNameAttribute_(const char* name) {
    if(!name) return;
    xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"name", "%s", name);
  }
  
  void   writeSizeAttribute_(unsigned int sz) {
    xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"size", "%u", sz);
  }
};







///////////////////////// INLINED CODE /////////////////////////
#define assert_reader_ok if(!reader_) throw XMLException( BuildExcMsg("no file opened") )

int xmlistream::nextOpeningTag_(xmlChar const* tagname) const throw(XMLException)
{
  assert_reader_ok;
  
  int ret, type;
  while(1) {
    ret = xmlTextReaderRead(reader_);
    if(ret==0) 
      return ret;
    if(ret<0)
      throw XMLException("unable to find next opening element");
    type = xmlTextReaderNodeType(reader_);
    if(type==XML_READER_TYPE_ELEMENT)
      break;
  }
  if(!tagname)
    return 1;
  
  xmlChar const* cur_tagname=xmlTextReaderConstName(reader_);
  if( xmlStrcmp(tagname,cur_tagname)==0 )
    return 1;
  throw XMLException("unable to find next opening element");
}


int xmlistream::nextClosingTag_(xmlChar const* tagname) const throw(XMLException)
{
  assert_reader_ok;
  
  int ret, type;
  while(1) {
    ret = xmlTextReaderRead(reader_);
    if(ret==0)
      return ret;
    if(ret<0)
      throw XMLException( BuildExcMsg("unable to find next closing tag") );
    type = xmlTextReaderNodeType(reader_);
    if(type==XML_READER_TYPE_END_ELEMENT)
      break;
  }
  if(!tagname) 
    return 1;
  
  xmlChar const* cur_tagname=xmlTextReaderConstName(reader_);
  if( xmlStrcmp(tagname, cur_tagname) == 0 )
    return 1;
  throw XMLException( BuildExcMsg("unable to find next closing tag") ); 
}


int xmlistream::nextTextElement_() const throw(XMLException)
{
  assert_reader_ok;
  
  int ret, type;
  while(1) {
    ret = xmlTextReaderRead(reader_);
    if(ret==0)
      return ret;
    if(ret<0) 
      throw XMLException( BuildExcMsg("unable to find next text element") );
    type = xmlTextReaderNodeType(reader_);
    if(type==XML_READER_TYPE_TEXT) 
      return 1; // should also add CDATA, maybe ?
  }
  throw XMLException( BuildExcMsg("unable to find next text element") );
}


int xmlistream::nextElement_(xmlChar const* tagname) const throw(XMLException)
{
  assert_reader_ok;
  
  int ret, return_value, type;
  while(1) {
    ret = xmlTextReaderRead(reader_);
    if(ret==0) 
      return ret;
    if(ret<0)
      throw XMLException("error while reading the next element");
    
    type = xmlTextReaderNodeType(reader_);
    if(type==XML_READER_TYPE_ELEMENT) {
      return_value=1;
      break;
    }
    if(type==XML_READER_TYPE_END_ELEMENT) {
      return_value=2;
      break;
    }
  }
  if(!tagname)
    return return_value;
  
  xmlChar const* cur_tagname=xmlTextReaderConstName(reader_);
  if( xmlStrcmp(tagname,cur_tagname)==0 )
    return return_value;
  throw XMLException("error while reading the next element");
}


void xmlistream::start_object(std::string& name, unsigned int& version) const
{
  assert_reader_ok;
  
  int ret;
  ret = nextOpeningTag_();
  readElementName_(name);
  readElementVersion_(version);
}


void xmlistream::end_object() const
{
  assert_reader_ok;
  nextClosingTag_();
}


void xmlistream::start_raw_pointer() const
{
  assert_reader_ok;
  nextOpeningTag_((xmlChar*)"pointer");
}


void xmlistream::end_raw_pointer() const
{
  assert_reader_ok;
  nextClosingTag_((xmlChar*)"pointer");
}


void xmlistream::start_reference() const
{
  assert_reader_ok;
  nextOpeningTag_((xmlChar*)"ref");
}


void xmlistream::end_reference() const
{
  assert_reader_ok;
  nextClosingTag_((xmlChar*)"ref");
}


void xmlistream::start_collection(unsigned int& size) const
{
  assert_reader_ok;
  
  nextOpeningTag_((xmlChar*)"collection");
  xmlChar* sz_val = xmlTextReaderGetAttribute(reader_, (xmlChar*)"size");
  size = atoi((char*)sz_val);
  xmlFree(sz_val);
}


void xmlistream::end_collection() const
{
  assert_reader_ok;
  nextClosingTag_((xmlChar*)"collection");
}




#define simple_type_read_def(type_, typname_, cvfunc_)  \
void xmlistream::read(type_& v) const                   \
{                                                       \
  assert_reader_ok;                                     \
  nextOpeningTag_();                                    \
  xmlChar* typ = xmlTextReaderGetAttribute(reader_, (xmlChar*)"t"); \
  if(typ==0) {                                          \
    std::string message = "empty value";                \
    return;                                             \
  }                                                     \
  if(xmlStrcmp((xmlChar*)#typname_,typ)!=0) {                       \
    std::string message = "gloups!";                    \
    return;                                             \
  }                                                     \
  xmlChar* value = xmlTextReaderGetAttribute(reader_, (xmlChar*)"v"); \
  if(value==0) {                                        \
    std::string message = "empty value";                \
    return;                                             \
  }                                                     \
  v = (type_)cvfunc_((const char*)value);               \
  xmlFree(typ);                                         \
  xmlFree(value);                                       \
}                                                       \

simple_type_read_def(int1,   i1,  atoi)
simple_type_read_def(uint1,  u1,  atoi)
simple_type_read_def(int2,   i2,  atoi)
simple_type_read_def(uint2,  u2,  atoi)
simple_type_read_def(int4,   i4,  atoi)
simple_type_read_def(uint4,  u4,  atoi)
simple_type_read_def(int8,   i8,  atol)
simple_type_read_def(uint8,  u8,  atol)
simple_type_read_def(float4, f4,  atof)
simple_type_read_def(float8, f8,  atof)
simple_type_read_def(std::string, string,)

#undef simple_type_read_def
#undef assert_reader_ok


////////////////////////////////////////////////////////////////////////////
#define assert_writer_ok if(!writer_) throw XMLException( BuildExcMsg("no file opened") )

void xmlostream::start_object(const std::string& name, unsigned int version)
{
  assert_writer_ok;
  //  xmlTextWriterWriteFormatString(writer_, "\n");
  xmlTextWriterStartElement(writer_, (xmlChar*)"object");
  xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"xsi:type",
				    "%s", (xmlChar*)name.c_str());
  xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"version",
				    "%d", version);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


void xmlostream::end_object()
{
  assert_writer_ok;
  xmlTextWriterEndElement(writer_);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


void xmlostream::start_raw_pointer(const char* name)
{
  assert_writer_ok;
  //  xmlTextWriterWriteFormatString(writer_, "\n");
  xmlTextWriterStartElement(writer_, (xmlChar*)"pointer");
  if(name) writeNameAttribute_(name);
}


void xmlostream::end_raw_pointer()
{
  assert_writer_ok;
  xmlTextWriterEndElement(writer_);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


void xmlostream::start_reference(const char* name)
{
  assert_writer_ok;
  //  xmlTextWriterWriteFormatString(writer_, "\n");
  xmlTextWriterStartElement(writer_, (xmlChar*)"ref");
  if(name) writeNameAttribute_(name);
}


void xmlostream::end_reference()
{
  assert_writer_ok;
  xmlTextWriterEndElement(writer_);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


void xmlostream::start_collection(unsigned int sz, const char* name)
{
  assert_writer_ok;
  //  xmlTextWriterWriteFormatString(writer_, "\n");
  xmlTextWriterStartElement(writer_, (xmlChar*)"collection");
  writeSizeAttribute_(sz);
  if(name) writeNameAttribute_(name);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


void xmlostream::end_collection()
{
  assert_writer_ok;
  xmlTextWriterEndElement(writer_);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


#define simple_type_write_def_OLD(type, tagname, format)       \
 void xmlostream::write(type v, const char* name)          \
 {                                                         \
   if(!writer_) {                                          \
     std::string message = "no file opened";               \
     return;                                               \
   }                                                       \
   xmlTextWriterStartElement(writer_, (xmlChar*)#tagname); \
   if(name)  writeNameAttribute_(name);                    \
   xmlTextWriterWriteFormatString(writer_, #format, v);    \
   xmlTextWriterEndElement(writer_);                       \
   xmlTextWriterWriteFormatString(writer_, "\n");          \
 }\


#define simple_type_write_def(type, typname, format)                         \
void xmlostream::write(type v, const char* name)                             \
{                                                                            \
   assert_writer_ok;                                                         \
   if(name==0) {                                                             \
     xmlTextWriterStartElement(writer_, (xmlChar*)"object");                 \
   }                                                                         \
   else {                                                                    \
     xmlTextWriterStartElement(writer_, (xmlChar*)name);                     \
   }                                                                         \
   xmlTextWriterWriteAttribute(writer_, (xmlChar*)"t", (xmlChar*)#typname);  \
   xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"v",#format,v);      \
   xmlTextWriterEndElement(writer_);                                         \
   xmlTextWriterWriteFormatString(writer_, "\n");                            \
}                                                                            \

simple_type_write_def(int1,   i1, %d)
simple_type_write_def(uint1,  u1, %u)
simple_type_write_def(int2,   i2, %d)
simple_type_write_def(uint2,  u2, %u)
simple_type_write_def(int4,   i4, %d)
simple_type_write_def(uint4,  u4, %u)
simple_type_write_def(int8,   i8, %ld)
simple_type_write_def(uint8,  u8, %lu)
simple_type_write_def(float4, f4, %.6E)
simple_type_write_def(float8, f8, %.12E)

void xmlostream::write(std::string const& v, const char* name)
{
  assert_writer_ok;
  if(name==0) {
    xmlTextWriterStartElement(writer_, (xmlChar*)"object");
  }
  else {
    xmlTextWriterStartElement(writer_, (xmlChar*)name);
  }
  xmlTextWriterWriteAttribute(writer_, (xmlChar*)"t", (xmlChar*)"string");
  xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"v","%s",v.c_str());
  xmlTextWriterEndElement(writer_);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


#undef simple_type_write_def
#undef assert_writer_ok


#endif

















#ifdef GARBAGE

#define simple_type_read_def_OLD(type_, tagname_, cvfunc_) \
  void xmlistream::read(type_& v) const \
  { \
    if(!reader_) {  \
      std::string message = "no file opened"; \
      return; \
    } \
    nextOpeningTag_((xmlChar*)#tagname_); \
    nextTextElement_(); \
    xmlChar const* value = xmlTextReaderConstValue(reader_); \
    if(value==0) { \
      std::string message = "empty value"; \
      return; \
    } \
    v = (type_)cvfunc_((const char*)value); \
    nextClosingTag_(); \
  } \

#endif

