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

#include "xmlschema_dict_writer.h"
#include "swig_dict_reader.h"



static struct option dictgen2_options[] = {
  {"list",             0, 0, 'l'},
  {"list-members",     0, 0, 'L'},
  {"all",              0, 0, 'a'},
  {"output",           1, 0, 'o'},
  {"xmlschema-output", 1, 0, 'x'},
  {"swig-output",      1, 0, 's'},
  {"verbose",          0, 0, 'v'},
  {"help",             0, 0, 'h'},
  {0, 0, 0, 0}
};



void usage()
{
  std::cerr << "usage: dictgen2 [OPTIONS] header_file                          " << std::endl;
  std::cerr << " where OPTIONS can be:                                         " << std::endl;
  std::cerr << "   --list,-l             list the persistent struct and classes" << std::endl;
  std::cerr << "   --list-members,-L     list the contents of the persistent objects" << std::endl;
  std::cerr << "   --all,-a              list also the non-persistent classes  " << std::endl;
  std::cerr << "   --output,-o           output file name                      " << std::endl;
  std::cerr << "   --swig-output,-s      swig original output file name        " << std::endl;
  std::cerr << "   --xmlschema-output,-x XMLSchema output file name            " << std::endl;
  std::cerr << "   --help,-h             print this message                    " << std::endl;
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
  int opt_idx=0;
  while( (c=getopt_long(argc, argv, "hvlLao:x:s:", 
			dictgen2_options, &opt_idx)) != -1 ) 
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
    case 'h':
      usage();
      break;
    default:
      usage();
    }
  if(optind>=argc) usage();
  headerName=argv[optind];

  std::string cmd;
  std::string tmp_name;  
  std::stringstream sstrm;
  sstrm << getpid();
  tmp_name = "dictgen_" + sstrm.str() + (std::string)".xml";  
  //  std::string cmd = "swig -c++ -xmllite -xml -module persist -o " 
  //    + tmp_name + " " + headerName;
  cmd = "swig -c++ -xml -w401 -module persist -o " 
    + tmp_name + " " + headerName;
  if(verbose)
    std::cout << cmd << std::endl;
  system(cmd.c_str());
  
  swig_dict_reader dr(tmp_name);
  xmlschema_dict_writer xsdw(xmlSchemaOutputFileName);
  codegen cg(headerName);
  if(generate_persisters)
    cg.openSourceFiles(outputFileName);
  
  dict d;
  int i,sz=dr.size();
  for(i=0;i<sz;i++) {
    dr.read(d);
    if(list && !list_members) 
      d.print(0);
    if(list_members) 
      d.print(1);
    if(generate_xmlschema)
      xsdw.write(d);
    if(generate_persisters)
      cg.generatePersister(d);
  }
  
  if(tmp_name!="" && !generate_swig_stuff) 
    remove(tmp_name.c_str());
  
  if(generate_swig_stuff)
    rename(tmp_name.c_str(), swigOutputFileName.c_str());
  
  exit(0);  

  //  xmlTextReaderPtr reader=0;

  
  
  // OK, now, we should have a xmlTextReader pointing on a SWIG XML file
  //  dict myDict;
  //  codegen cg(headerName);
  //  if(outputFileName!="")
  //    cg.openSourceFiles(outputFileName);
  
  //  dict_output<xmlostream> dout("test_swig.dict",0);
  //  while( /*myDict.readFromSwigXMLStream(reader) ==*/1 ) {
  //    if(list)
  //      myDict.print(verbose);
  //    if(generate_persisters)
  //      cg.generatePersister(myDict);
  //    if(generate_xmlschema) {
  //      std::string filename = xmlSchemaOutputFileName.substr(0, xmlSchemaOutputFileName.find(".xml"));
  //      filename = filename.substr(0, filename.find(".xsd"));
  //      filename = filename + "_" + myDict.name() + ".xsd";
  //      std::cout << " generate_schema ! filename=" << filename << std::endl;
  //    }
  //  }
  
  
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

