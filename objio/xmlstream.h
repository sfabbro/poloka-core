// -*- C++ -*-
// 
// $Id: xmlstream.h,v 1.5 2004/03/02 12:40:16 nrl Exp $
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



class xmlstream {
public:
  xmlstream();
  xmlstream(std::string const& filename, int mode, int compression);
  ~xmlstream();
  
  void            open(std::string const&, int mode, int compression);
  void            close();
  
  inline void     read_start_object_tag(std::string& name, unsigned int& version) const;
  inline void     read_end_object_tag() const;
  
  inline void     read_start_raw_pointer_tag() const;
  inline void     read_end_raw_pointer_tag() const;
  
  inline void     read_start_collection_tag(unsigned int&) const;
  inline void     read_end_collection_tag() const;
  
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
  
  template<class T>
  inline void     read(T*&) const;
  
  inline void     skip() const;
  
  
  inline void     write_start_object_tag(std::string const& name, unsigned int version, void const* addr=0);
  inline void     write_end_object_tag();
  
  inline void     write_start_raw_pointer_tag(const char* name, void const* addr=0);
  inline void     write_end_raw_pointer_tag();
  
  inline void     write_start_collection_tag(unsigned int size, const char* name=0, void const* addr=0);
  inline void     write_end_collection_tag();
  
  inline void     write(int1 v, const char* name=0,   void const* addr=0);
  inline void     write(uint1 v, const char* name=0,  void const* addr=0);
  inline void     write(int2 v, const char* name=0,   void const* addr=0);
  inline void     write(uint2 v, const char* name=0,  void const* addr=0);
  inline void     write(int4 v, const char* name=0,   void const* addr=0);
  inline void     write(uint4 v, const char* name=0,  void const* addr=0);
  inline void     write(int8 v, const char* name=0,   void const* addr=0);
  inline void     write(uint8 v, const char* name=0,  void const* addr=0);
  inline void     write(float4 v, const char* name=0, void const* addr=0);
  inline void     write(float8 v, const char* name=0, void const* addr=0);
  inline void     write(const std::string& v, const char* name=0, void const* addr=0);
  
  template<class T>
  inline void     write(T const*, const char* name=0);
  
private:
  xmlTextReaderPtr reader_;  
  xmlTextWriterPtr writer_;
  
  inline int nextOpeningTag_(xmlChar const* tagname=0) const throw(XMLException);
  inline int nextClosingTag_(xmlChar const* tagname=0) const throw(XMLException);
};



///////////////////////// INLINE CODE /////////////////////////
#define assert_reader_ok if(!reader_) throw XMLException( BuildExcMsg("no file opened") )
#define assert_writer_ok if(!writer_) throw XMLException( BuildExcMsg("no file opened") )



int xmlstream::nextOpeningTag_(xmlChar const* tagname) const throw(XMLException)
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


int xmlstream::nextClosingTag_(xmlChar const* tagname) const throw(XMLException)
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



void xmlstream::read_start_object_tag(std::string& name, unsigned int& version) const
{
  assert_reader_ok;
  
  nextOpeningTag_();
  xmlChar const* nm = xmlTextReaderConstName(reader_);
  name = (char*)nm;
  version=0;
  xmlChar* ver_str = xmlTextReaderGetAttribute(reader_, (xmlChar*)"version");
  if(ver_str==0) return; // should throw an exc. here
  version = (unsigned int)atoi((char*)ver_str);
}


void xmlstream::read_end_object_tag() const
{
  assert_reader_ok;
  nextClosingTag_();
}


void xmlstream::read_start_raw_pointer_tag() const
{
  assert_reader_ok;
  nextOpeningTag_((xmlChar*)"pointer");
}


void xmlstream::read_end_raw_pointer_tag() const
{
  assert_reader_ok;
  nextClosingTag_((xmlChar*)"pointer");
}


void xmlstream::read_start_collection_tag(unsigned int& size) const
{
  assert_reader_ok;
  
  nextOpeningTag_((xmlChar*)"collection");
  xmlChar* sz_val = xmlTextReaderGetAttribute(reader_, (xmlChar*)"size");
  size = atoi((char*)sz_val);
  xmlFree(sz_val);
}


void xmlstream::read_end_collection_tag() const
{
  assert_reader_ok;
  nextClosingTag_((xmlChar*)"collection");
}


#define simple_type_read_def(type_, typname_, cvfunc_)  \
void xmlstream::read(type_& v) const                    \
{                                                       \
  assert_reader_ok;                                     \
  nextOpeningTag_();                                    \
  xmlChar* typ = xmlTextReaderGetAttribute(reader_, (xmlChar*)"xsi:type"); \
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

void xmlstream::skip() const
{
  assert_reader_ok;
  nextOpeningTag_();
  nextClosingTag_();
}



void xmlstream::write_start_object_tag(const std::string& name, unsigned int version, void const* addr)
{
  assert_writer_ok;
  
  xmlTextWriterStartElement(writer_, (xmlChar*)"object");
  xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"xsi:type",
				    "%s", (xmlChar*)name.c_str());
  xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"version",
				    "%d", version);
  if(addr) xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"addr",
					     "%lu", addr);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


void xmlstream::write_end_object_tag()
{
  assert_writer_ok;
  xmlTextWriterEndElement(writer_);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


void xmlstream::write_start_raw_pointer_tag(const char* name, void const* addr)
{
  assert_writer_ok;
  xmlTextWriterStartElement(writer_, (xmlChar*)"pointer");
  if(name) xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"name", "%s", name);
  if(addr) xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"addr",
					     "%lu", addr);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


void xmlstream::write_end_raw_pointer_tag()
{
  assert_writer_ok;
  xmlTextWriterEndElement(writer_);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


void xmlstream::write_start_collection_tag(unsigned int sz, const char* name, void const* addr)
{
  assert_writer_ok;
  xmlTextWriterStartElement(writer_, (xmlChar*)"collection");
  xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"size", "%u", sz);
  if(name) xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"name", "%s", name);
  if(addr) xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"addr",
					     "%lu", addr);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


void xmlstream::write_end_collection_tag()
{
  assert_writer_ok;
  xmlTextWriterEndElement(writer_);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


#define simple_type_write_def(type, typname, format)                         \
void xmlstream::write(type v, const char* name, void const* addr)            \
{                                                                            \
   assert_writer_ok;                                                         \
   if(name==0) {                                                             \
     xmlTextWriterStartElement(writer_, (xmlChar*)"object");                 \
   }                                                                         \
   else {                                                                    \
     xmlTextWriterStartElement(writer_, (xmlChar*)name);                     \
   }                                                                         \
   xmlTextWriterWriteAttribute(writer_, (xmlChar*)"xsi:type", (xmlChar*)#typname);  \
   xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"v",#format,v);      \
   if(addr) xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"addr",     \
                                              "%lu", addr);                  \
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

void xmlstream::write(std::string const& v, const char* name, void const* addr)
{
  assert_writer_ok;
  if(name==0) {
    xmlTextWriterStartElement(writer_, (xmlChar*)"object");
  }
  else {
    xmlTextWriterStartElement(writer_, (xmlChar*)name);
  }
  xmlTextWriterWriteAttribute(writer_, (xmlChar*)"xsi:type", (xmlChar*)"string");
  xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"v","%s",v.c_str());
  if(addr) xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"addr",
					     "%lu", addr);
  xmlTextWriterEndElement(writer_);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


template<class T>
inline void xmlstream::write(T const* t, const char* name)
{
  assert_writer_ok;
  
  if(name==0)
    xmlTextWriterStartElement(writer_, (xmlChar*)"address");
  else {
    xmlTextWriterStartElement(writer_, (xmlChar*)name);
  }
  xmlTextWriterWriteFormatAttribute(writer_, (xmlChar*)"v","%lu",(void*)t);
  xmlTextWriterEndElement(writer_);
  xmlTextWriterWriteFormatString(writer_, "\n");
}


//template<class T>
//void xmlstream::write(T const* t)
//{
//  assert_writer_ok;
//  if(name==0) {
//    
//  }
//}

#undef simple_type_read_def
#undef assert_reader_ok
#undef simple_type_write_def
#undef assert_writer_ok


#endif

