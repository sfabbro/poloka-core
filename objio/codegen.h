// -*- C++ -*-
// 
// file codegen.h
// 
// 
#ifndef CODEGEN_H
#define CODEGEN_H

#include <iostream>
#include <fstream>

#include "dict.h"


class codegen {
public:
  codegen(std::string const& header_name) 
    : source_header_name_(header_name) {}
  ~codegen() { closeSourceFiles(); }
  
  void  openSourceFiles(std::string const& base);
  void  openClassHeaderFile(std::string const& className);
  void  generatePersister(dict const&);
  void  closeSourceFiles();
  void  closeClassHeaderFile();
  
private:
  
  void  classPersisterDecl_(dict const&);
  void  classPersisterCode_(dict const&);
  
  std::string   source_header_name_;
  std::ofstream ofs_h_;
  std::ofstream ofs_cc_;
  std::ofstream class_ofs_h_;
  
  std::vector<std::string> dependencies_;
};



#endif


