// -*- C++ -*-
// 
// file dictgen2.cc
// 
//
#define _GNU_SOURCE
#include <getopt.h>
 
#include <sys/types.h>
#include <unistd.h>

#include <iostream>
#include <sstream>
#include <string>

#include "dict.h"
#include "xmlstream.h"
#include "dictio.h"
#include "codegen.h"


static struct option dictgen2_options[] = {
  {"list",             0, 0, 'l'},
  {"list-members",     0, 0, 'L'},
  {"all",              0, 0, 'a'},
  {"output",           1, 0, 'o'},
  {"xmlschema-output", 1, 0, 'x'},
  {"swig-output",      1, 0, 's'},
  {"verbose",          0, 0, 'v'},
  {0, 0, 0, 0}
};



void usage()
{
  std::cerr << "usage: dictgen2 [OPTIONS] header_file                          " << std::endl;
  std::cerr << " where OPTIONS can be:                                         " << std::endl;
  std::cerr << "   --list,-l             list the persistent struct and classes" << std::endl;
  std::cerr << "   --list-members,-L     list the contents of the persistent   " << std::endl
	    << "                            as well as their contents          " << std::endl;
  std::cerr << "   --all,-a              list also the non-persistent classes  " << std::endl;
  std::cerr << "   --output,-o           output file names                     " << std::endl;
  std::cerr << "   --swig-output,-s      swig original output                  " << std::endl;
  std::cerr << "   --xmlschema-output,-x XMLSchema output file name            " << std::endl;
  exit(-1);
}



int 
main(int argc, char** argv)
{
  std::string xmlInputFileName;
  std::string xmlOutputFileName;
  
  int list=0;
  int verbose=0;
  int list_members=0;
  int generate_persisters=0;
  int generate_xmlschema=0;
  int generate_swig_stuff=0;
  
  std::string headerName;
  std::string outputFileName="";
  std::string swigOutputFileName="";
  std::string xmlSchemaOutputFileName="";
  
  
  char c;
  while( (c=getopt(argc, argv, "vlLao:x:s:")) != -1 ) 
    switch(c) {
    case 'v':
      verbose=1;
      break;
    case 'l':
      list=1;
      break;
    case 'L':
      list=1;
      list_members=1;
      break;
    case 'o':
      outputFileName=optarg;
      generate_persisters=1;
      break;
    case 'x':
      generate_xmlschema=1;
      xmlSchemaOutputFileName=optarg;
      break;
    case 's':
      generate_swig_stuff=1;
      swigOutputFileName=optarg;
      break;
    default:
      usage();
    }
  if(optind>=argc) usage();
  headerName=argv[optind];
  exit(0);
  
  
  xmlTextReaderPtr reader=0;
  std::string tmp_name;  
  
  // if no XML input file specified, 
  // we have to generate it from the header...
  if(xmlInputFileName=="") {
    if(optind>=argc) usage();
    headerName=argv[optind];
    
    // here, we try to run swig on the header file...
    if(xmlOutputFileName=="") {
      std::stringstream sstrm;
      sstrm << getpid();
      tmp_name = "dictgen_" + sstrm.str() + (std::string)".xml";
      xmlOutputFileName = tmp_name;
    }
    std::string cmd = "swig -c++ -xmllite -xml -module persist -o " 
      + tmp_name + " " + headerName;
    std::cout << cmd << std::endl;
    system(cmd.c_str());
    reader = xmlNewTextReaderFilename(tmp_name.c_str());
    if(reader==0) {
      std::cout << " ERROR no file opened !" << std::endl;
    }
  }
  // a XML input file was provided. Just have to open it.
  else {
    reader = xmlNewTextReaderFilename(xmlInputFileName.c_str());
  }
  
  //  if(outputFileName=="") {
  //    if(headerName!="") 
  //      outputFileName = headerName.substr(0,headerName.find("."));
  //    else
  //      outputFileName = "dictgen";
  //  }
  //  else {
  //    outputFileName = outputFileName.substr(0,outputFileName.find(".cc"));
  //    outputFileName = outputFileName.substr(0,outputFileName.find(".h"));    
  //  }
  
  // OK, now, we should have a xmlTextReader pointing on a SWIG XML file
  dict myDict;
  codegen cg(headerName);
  if(outputFileName!="")
    cg.openSourceFiles(outputFileName);
  
  //  dict_output<xmlostream> dout("test_swig.dict",0);
  while( /*myDict.readFromSwigXMLStream(reader) ==*/1 ) {
    if(list)
      myDict.print(verbose);
    if(generate_persisters)
      cg.generatePersister(myDict);
    if(generate_xmlschema) {
      std::string filename = xmlSchemaOutputFileName.substr(0, xmlSchemaOutputFileName.find(".xml"));
      filename = filename.substr(0, filename.find(".xsd"));
      filename = filename + "_" + myDict.name() + ".xsd";
      std::cout << " generate_schema ! filename=" << filename << std::endl;
    }
  }
  
  if(tmp_name!="") 
    remove(tmp_name.c_str());
  
  //  dout.close();
}






//  xmlostream ostrm_tst;
//  std::list<std::string> l;
//  l.push_back("toto");
//  l.push_back("tutu");
//  l.push_back("titi");
//  ostrm_tst.open("toto.xsd");
//  ostrm_tst.start_dict("toto", "struct", 10, 2, l);
//  ostrm_tst.write_element("x_", "f4");
//  ostrm_tst.write_element("y_", "f4");
//  ostrm_tst.write_element("z_", "f4");
//  ostrm_tst.write_element("id_", "u4");
//  ostrm_tst.end_dict();
//  ostrm_tst.close();
//  exit(0);

