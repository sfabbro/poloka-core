// -*- C++ -*-
// 
// file codegen.cc
// 
// 
#include <iostream>
#include <string.h>

#include "codegen.h"


void codegen::openSourceFiles(std::string const& base)
{
  std::string target_header_name = base + ".h";
  ofs_h_.open(target_header_name.c_str());
  ofs_h_ << "// -*- C++ -*-" << std::endl
	 << "#ifndef __" << (target_header_name.substr(0,target_header_name.find("."))) << "__" << std::endl
	 << "#define __" << (target_header_name.substr(0,target_header_name.find("."))) << "__" << std::endl
	 << std::endl
	 << "#include <string>" << std::endl
	 << "#include \"persister.h\"" << std::endl
	 << "#include \"objio.h\"" << std::endl
	 << "#include \"" << source_header_name_ << "\"" << std::endl
	 << std::endl;
  
  std::string target_code_name = base + ".cc";
  ofs_cc_.open(target_code_name.c_str());
  ofs_cc_ << "// -*- C++ -*-" << std::endl
	  << "#include \"" << target_header_name.c_str() << "\"" << std::endl
	  << std::endl;
}


void codegen::closeSourceFiles()
{
  if(ofs_h_.is_open()) {
    ofs_h_ << std::endl
	   << "#endif" << std::endl;
    ofs_h_.close();
  }
  if(ofs_cc_.is_open())
    ofs_cc_.close();
}


void codegen::generatePersister(dict const& d)
{
  classPersisterDecl_(d);
  classPersisterCode_(d);
}


void codegen::classPersisterDecl_(dict const& dict_)
{
  //  if(!dict_.isPersistent()) return;
  
  ofs_h_ << "template<>" << std::endl
	 << "class persister<" << dict_.name() << "> : public handle<" 
	 << dict_.name() << "> {" << std::endl
	 << "public:" << std::endl
    // constructors & destructor
	 << "persister<" << dict_.name() << ">() : handle<" 
	 << dict_.name() << ">() {}" << std::endl
	 << "persister<" << dict_.name() << ">(" 
	 << dict_.name() << "& obj) : handle<" 
	 << dict_.name() << ">(&obj) {}" << std::endl
	 << "persister<" << dict_.name() << ">(" 
	 << dict_.name() << " const & obj) : handle<" 
	 << dict_.name() << ">(const_cast<" 
	 << dict_.name() << "*>(&obj)) {}// I know, that's UGLY." << std::endl
	 << "~persister<" << dict_.name() 
	 << ">() {}" << std::endl << std::endl
    
    // interface methods
	 << "unsigned int version() const { return " << dict_.version() << "; }" << std::endl
	 << "std::string  name()    const { return (std::string)\"" << dict_.name() << "\"; }" << std::endl
	 << "unsigned int size()    const { return " << dict_.size() << ";}" << std::endl
	 << "std::string  name(unsigned int i) const { return (std::string)memberNames_[i]; }" << std::endl
	 << "std::string  type(unsigned int i) const { return (std::string)memberTypes_[i]; }" << std::endl
	 << std::endl
	 << "private:" << std::endl
	 << "template<class IOS>" << std::endl
	 << "void write_members(obj_output<IOS>& oo) const {" << std::endl;
  
  
  unsigned int i;
  for(i=0;i<dict_.size();i++) 
    if(dict_.member(i).isPersistent)
      ofs_h_ << "oo.write(obj_->"  << dict_.member(i).name  
	     << ", \""  << dict_.member(i).name  
	     << "\"" << ");" << std::endl;
  
  ofs_h_ << "}" << std::endl
	 << std::endl
	 << "template<class IOS>"
	 << "void read_members(obj_input<IOS> const& oi) {" << std::endl;
  
  for(i=0;i<dict_.size();i++)
    if(dict_.member(i).isPersistent)
      ofs_h_ << "oi.read(obj_->"  << dict_.member(i).name  
	     << ");" << std::endl;
  
  ofs_h_ << "}" << std::endl
	 << std::endl
	 << "template<class T> friend class obj_input;" << std::endl
	 << "template<class T> friend class obj_output;" << std::endl
	 << std::endl
	 << "static const char* memberNames_[];" << std::endl
	 << "static const char* memberTypes_[];" << std::endl
	 << "};" << std::endl
	 << std::endl
	 << "template<class IOS>" << std::endl
	 << "obj_output<IOS>& operator<<(obj_output<IOS>& oo, " << dict_.name() << " const& p)" << std::endl
	 << "{" << std::endl
	 << "persister<" << dict_.name() << "> pp(p);" << std::endl
	 << "oo.write(pp);" << std::endl
	 << "return oo;" << std::endl
	 << "}" << std::endl
	 << std::endl
	 << "template<class IOS>" << std::endl
	 << "obj_input<IOS> const& operator>>(obj_input<IOS> const& oi, " << dict_.name() << "& p)" << std::endl
	 << "{" << std::endl
	 << "persister<" << dict_.name() << "> pp(p);" << std::endl
	 << "oi.read(pp);" << std::endl
	 << "return oi;" << std::endl
	 << "};" << std::endl;
}




void codegen::classPersisterCode_(dict const& dict_)
{
  if(!dict_.isPersistent()) return;
  
  ofs_cc_ << "const char* persister<" << dict_.name() << ">::memberNames_[] = {";
  if(dict_.size()>0) {
    unsigned int i;
    ofs_cc_ << "\""  << dict_.member(0).name  << "\"";
    for(i=1;i<dict_.size();i++)
      ofs_cc_ << ", \""  << dict_.member(i).name  << "\"";
  }
  ofs_cc_ << "};" << std::endl;
  
  ofs_cc_ << "const char* persister<" << dict_.name() << ">::memberTypes_[] = {";
  if(dict_.size()>0) {
    unsigned int i;
    ofs_cc_ << "\""  << dict_.member(0).type  << "\"";
    for(i=1;i<dict_.size();i++)
      ofs_cc_ << ", \""  << dict_.member(i).type  << "\"";
  }
  ofs_cc_ << "};" << std::endl;
}


