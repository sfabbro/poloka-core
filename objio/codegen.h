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
#include "cppclass.h"


class codegen {
public:
  codegen(std::string const& header_name, std::string const& config_file="")
    : counter_(0), source_header_name_(header_name), config_file_name_(config_file) {}
  ~codegen() { closeOutputFiles(); }
  
  void  openOutputFiles(std::string const& base);
  void  closeOutputFiles();
  
  void  openClassHeaderFile(std::string const& className);
  void  closeClassHeaderFile();
  
  void  generatePersister(dict const&);
  void  generatePersister(CppClass const&);
  
  static std::string cleanTypeName(std::string const&);
  static std::string cleanTemplateArgList(CppClass const&);
  
private:

  void  classPersisterDecl_(CppClass const&);
  void  classPersisterCode_(CppClass const&);
  
  void  classPersisterDecl_(dict const&);
  void  classPersisterCode_(dict const&);
  
  bool checkTemplateInstantiation_(CppClass const&,
				   CppClass& checkedClass,
				   std::vector<CppClass>& classVec);
  
  int findNextCket(std::string const& buffer,int after);
  int findNextComaOutsideBracket(std::string const& buffer,int after);
  int findRecursivelyTemplateInstantiation_(std::string & sbuff,
					    std::string const& className,
					    std::vector<templateInstantiation>& tvec);
  int readOneTemplateInstantiation_(std::string const& classname,
				    std::string & content,
				    std::vector<templateInstantiation>& tvec);
  void  checkTemplateInstantiation_(std::string const& className_,
				    std::vector<templateInstantiation>&);
  
  int  counter_;
  std::string   source_header_name_;
  std::string   config_file_name_;
  std::ofstream ofs_h_;
  std::ofstream ofs_cc_;
  std::ofstream class_ofs_h_;
};



#endif

