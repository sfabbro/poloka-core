// -*- C++ -*-
// 
// file xmlstream.cc
// 
// 
#include <assert.h>

#include <list>

#include <config.h>

#include "xmlstream.h"
#include "objio_defs.h"


xmlstream::xmlstream()
  : reader_(0), writer_(0)
{
}


xmlstream::xmlstream(std::string const& filename, int mode, int compression)
  : reader_(0), writer_(0)
{
  open(filename, mode, compression);
}


xmlstream::~xmlstream()
{
  close();
}



void  xmlstream::open(std::string const& filename, int mode, int compression)
{
  close();
  
  if(mode==OBJIO_READ) {
    reader_ = xmlNewTextReaderFilename(filename.c_str());
    nextOpeningTag_(); //doc
    return;
  }
  if(mode==OBJIO_WRITE) {
    writer_ = xmlNewTextWriterFilename(filename.c_str(), compression);
    //    xmlCharEncodingHandlerPtr enc = xmlFindCharEncodingHandler("UTF-8");
    //    xmlOutputBufferPtr ptr = xmlAllocOutputBuffer(enc);
    //    writer_ = xmlNewTextWriter(ptr);
    xmlTextWriterStartDocument(writer_, "1.0", "UTF-8", "no");
    
    xmlTextWriterWriteFormatString(writer_, "\n\n\n");  
    
    xmlTextWriterWriteFormatComment(writer_, "\nwritten by %s version %s\n",
				    PACKAGE, VERSION);
    
    xmlTextWriterWriteFormatString(writer_, "\n\n\n");
    xmlTextWriterStartElement(writer_, (xmlChar*)"objects");
    xmlTextWriterWriteAttribute(writer_, (xmlChar*)"xmlns:xsi",
				(xmlChar*)"http://www.w3.org/2001/XMLSchema-instance");
    xmlTextWriterWriteAttribute(writer_, (xmlChar*)"xmlns:snls",
				(xmlChar*)"http://supernovae.in2p3.fr/SNLS");
    xmlTextWriterWriteFormatString(writer_, "\n");
  }
}


void  xmlstream::close()
{
  if(reader_) {
    xmlFreeTextReader(reader_);
    reader_=0;
  }
  if(writer_) {
    xmlTextWriterEndDocument(writer_);
    xmlFreeTextWriter(writer_);
    writer_=0;
  }
}

