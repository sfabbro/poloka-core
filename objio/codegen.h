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
  codegen(std::string const& header_name, std::string const& config_file="")
    : counter_(0), source_header_name_(header_name), config_file_name_(config_file) {}
  ~codegen() { closeOutputFiles(); }
  
  void  openOutputFiles(std::string const& base);
  void  closeOutputFiles();
  
  void  openClassHeaderFile(std::string const& className);
  void  closeClassHeaderFile();
  
  void  generatePersister(dict const&);
  
private:
  
  //  // this struct just tells us how to instantiate a template
  //  struct TemplateInst_ {
  //    std::vector<std::string> symbolicTypes_;
  //    std::vector<std::string> realTypes_;
  //    std::string fullName_;
  //    
  //    std::string    concatenateSymbolicTypes() const;
  //    std::string    concatenateRealTypes() const;
  //    
  //    void           copy(TemplateInst_ const& ti) {
  //      std::copy(ti.symbolicTypes_.begin(), ti.symbolicTypes_.end(), std::back_inserter(symbolicTypes_));
  //      std::copy(ti.realTypes_.begin(), ti.realTypes_.end(), std::back_inserter(realTypes_));
  //      fullName_ = ti.fullName_;
  //    }
  //    TemplateInst_& operator=(TemplateInst_ const& ti) { copy(ti); return *this; }
  //  };
  
  //  bool  parseCommandLine_(std::string const& buff, std::string const& rqst_className, TemplateInst_& ti);
 
  void  classPersisterDecl_(dict const&);
  void  classPersisterCode_(dict const&);

  
  int findNextCket(std::string const& buffer,int after);
  int findNextComaOutsideBracket(std::string const& buffer,int after);
  int findRecursivelyTemplateInstantiation_(std::string & sbuff,std::string const& className, std::vector<templateInstantiation>& tvec);
  int readOneTemplateInstantiation_(std::string const& classname,std::string & content, std::vector<templateInstantiation>& tvec);
  void  checkTemplateInstantiation_(std::string const& className_,std::vector<templateInstantiation>&);
  
  
  int  counter_;
  std::string   source_header_name_;
  std::string   config_file_name_;
  std::ofstream ofs_h_;
  std::ofstream ofs_cc_;
  std::ofstream class_ofs_h_;
  
  //  std::vector<std::string> dependencies_;
};



#endif

