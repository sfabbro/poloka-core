// -*- C++ -*-
// 
// file testxml.cc
// 
// 
#include <iostream>

#include <libxml/xmlreader.h>



int
main(int argc, char** argv)
{
  int ret, count=0;
  
  xmlTextReaderPtr reader = 
    xmlNewTextReaderFilename("test.xml");
  if(reader==0) {
    std::cout << "testxml: unable to open test.xml"
	      << std::endl;
    return -1;
  }
  
  while(1) {
    ret = xmlTextReaderRead(reader);
    if(ret==0) break;
    if(ret<0) {
      std::cout << "testxml: error while reading test.xml: count="
		<< count << " (ret=" << ret << ")" << std::endl;
    }
    count++;
  }
  
  std::cout << "testxml: " << count << " elements read." << std::endl;
}
