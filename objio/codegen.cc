// -*- C++ -*-
// 
// file codegen.cc
// 
// 
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include "codegen.h"

using namespace std;

void codegen::openOutputFiles(std::string const& base)
{
  std::string target_header_name = base + ".h";
  ofs_h_.open(target_header_name.c_str());
  ofs_h_ << "// -*- C++ -*-" << std::endl
	 << "#ifndef __" << (target_header_name.substr(0,target_header_name.find("."))) << "__H__" << std::endl
	 << "#define __" << (target_header_name.substr(0,target_header_name.find("."))) << "__H__" << std::endl
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
	  << "#include \"typemgr.h\"" << std::endl
	  << "#include \"xmlstream.h\"" << std::endl
	  << std::endl;
}


void codegen::closeOutputFiles()
{
  if(ofs_h_.is_open()) {
    ofs_h_ << std::endl
	   << "#endif" << std::endl;
    ofs_h_.close();
  }
  if(ofs_cc_.is_open())
    ofs_cc_.close();
}



void  codegen::openClassHeaderFile(std::string const& className)
{
  if(class_ofs_h_.is_open())
    class_ofs_h_.close();
  
  std::string class_header_name = className + "__persister.h";
  class_ofs_h_.open(class_header_name.c_str());
  
  class_ofs_h_ << "// -*- C++ -*-" << std::endl
	       << "#ifndef __" << className << "__PERSISTER__H__" << std::endl
	       << "#define __" << className << "__PERSISTER__H__" << std::endl
	       << std::endl
	       << "#include <string>" << std::endl
	       << "#include \"persister.h\"" << std::endl
	       << "#include \"objio.h\"" << std::endl
	       << "#include \"" << source_header_name_ << "\"" << std::endl
	       << std::endl;
}


void codegen::closeClassHeaderFile()
{
  if(class_ofs_h_.is_open()) {
    class_ofs_h_ << std::endl
	   << "#endif" << std::endl;
    class_ofs_h_.close();
  }
}


void codegen::generatePersister(dict const& d)
{
  std::vector<templateInstantiation> tvec;
  //checkConfigFile_(d.name(), tvec);
  checkTemplateInstantiation_(d.name(), tvec);
  
  //  std::cout << " checkConfigFile_... tvec.size()="
  //	    << tvec.size() << std::endl;
  
  //  int ii,sz=(int)tvec.size();
  //  for(ii=0;ii<sz;ii++)
  //    std::cout << tvec[ii].fullSymbolicName() 
  //	      << " --> " 
  //	      << tvec[ii].fullRealName()
  //	      << std::endl;
  
  openClassHeaderFile(d.name());
  if(!d.isTemplate()) {
    classPersisterDecl_(d);
    classPersisterCode_(d);
  }
  else {
    dict td = d;
    if(tvec.size()>0)
      td.update(tvec[0]);
    classPersisterDecl_(td);
    int i;
    for(i=0;i<(int)tvec.size();i++) {
      td = d;
      td.instantiate(tvec[i]);
      classPersisterCode_(td);
    }
  }
  closeClassHeaderFile();
}

// find next ',' index outside a bracket
// , <,<,,>> 
int codegen::findNextComaOutsideBracket(std::string const& buffer,int after) {
  int size = buffer.size();
  int pos_coma = buffer.find(',',after);
  if(pos_coma<0 || pos_coma > size)
    return -1;
  int pos_bra = buffer.find('<',after);
  if(pos_bra<0 || pos_bra > size) // no bracket, that's easy
    return pos_coma;
  if(pos_bra>pos_coma)
    return pos_coma; // we have something like ".. , ..<" so it's ok
  int pos_cket = findNextCket(buffer,pos_bra);
  if(pos_cket>pos_coma) { // we have something like " ...< , >" so we skip this coma
    return findNextComaOutsideBracket(buffer,pos_cket);
  }
  return pos_coma; // ... < > ,
}

// find in a string the '>' index matching with '<' at index after
// (this needs a function in case of < <> < <> > > ...) 
int codegen::findNextCket(std::string const& buffer,int after) {
  int nbrackets=1;
  int size = buffer.size();
  int index = after;
  while(index<size) {
    char c = buffer[index];
    if(c=='>') {
      nbrackets--;
    }
    if(nbrackets==0)
      return index;
    if(c=='<')
      nbrackets++;
    index++;
  }
  return -1;
}

// parse classname<content> which may be classname<c1,c2,c3<c4,c5>> ....
int codegen::readOneTemplateInstantiation_(std::string const& classname,std::string & content, std::vector<templateInstantiation>& tvec) {
  
  templateInstantiation ti;
  ti.templateClassName_ = classname;
  
  char symbolicstring[1];
  char symbolicchar = 'T';
  
  
  // split content with ","
  int pos = 0;
  while(true){
    pos = findNextComaOutsideBracket(content,pos);
    if(pos<0 || pos > content.size()) { 
      string oneclass = content;
      ti.realTypes_.push_back(content);
      sprintf(symbolicstring,"%c",symbolicchar++);
      ti.symbolicTypes_.push_back(symbolicstring); 
      tvec.push_back(ti);
      return 0;
    }
    string oneclass = content;
    oneclass.resize(pos);
    content.erase(0,pos+1);
    ti.realTypes_.push_back(oneclass);
    sprintf(symbolicstring,"%c",symbolicchar++);
    ti.symbolicTypes_.push_back(symbolicstring);
  }
  std::cout << " HELLO ?!" << std::endl;
  tvec.push_back(ti);
  return 0;
}


// find recusively templates in string sbuff, and save them
// return position of >

int codegen::findRecursivelyTemplateInstantiation_(std::string & sbuff, std::string const& className, std::vector<templateInstantiation>& tvec) {
  
  int size = sbuff.size();
  
  int newbegin = sbuff.find("<");
  if(newbegin<0 || newbegin>=size)
    return 0; // no bracket in sbuff
  int newend = findNextCket(sbuff,newbegin+1);
  
  // we have found matching < ..... > with indexes newbegin and new end
  // do stuff ... 
  
  // read classname
  string classnamestring = sbuff;
  classnamestring.resize(newbegin);
  stringstream stream(classnamestring.c_str());
  stream >> classnamestring;
  
  // check if classnamestring contains required className
  if(classnamestring==className) {
    // what is inside <...>
    string content = sbuff;
    content.erase(newend,content.size()-(newend));
    content.erase(0,newbegin+1);
    
    // read content
    readOneTemplateInstantiation_(className,content,tvec);
  }

  // remove beginning of sbuff already processed
  sbuff.erase(0,newend+1);
  
  // known try to find and other pair of matching  < ..... >
  findRecursivelyTemplateInstantiation_(sbuff,className,tvec);
  
  return 1;
 }

void codegen::checkTemplateInstantiation_(std::string const& className, std::vector<templateInstantiation>& tvec)
{
  //std::cout << std::endl << "codegen::checkTemplateInstantiation_ looking for " << className << " ..." << std::endl;
  tvec.clear();
  ifstream ifs;
  ifs.open(source_header_name_.c_str());
  if(!ifs.is_open()) {
      std::cout << "codegen::checkTemplateInstantiation_ ERROR, unable to open "
		<< source_header_name_ << std::endl;
      return;
  }
  
  // reading the source_header and look for lines containing key word "make_persister_for"
  char buff[2048];
  string sbuff;
  string keyword = "make_persister_for";
  
  // loop on lines
  while( ifs.good() ) {
    ifs.getline(buff, 2048);
    sbuff=buff;
    int pos = sbuff.find(keyword);
    if(pos<0 || pos>=sbuff.size())
      continue;
    int pos2 = sbuff.find(className);
    if(pos2<0 || pos2>=sbuff.size())
      continue;
    sbuff.erase(pos,keyword.size());
    
    pos = sbuff.find("//");
    if(pos>=0 && pos<sbuff.size())
      sbuff.erase(pos,2);
    
    findRecursivelyTemplateInstantiation_(sbuff,className,tvec);
  }
}


/* use checkTemplateInstantiation_ instead
   void codegen::checkConfigFile_(std::string const& className, std::vector<templateInstantiation>& tvec)
   {
   tvec.clear();
   ifstream ifs;
  
   // config file ? 
   if(config_file_name_!="") {
   ifs.open(config_file_name_.c_str());
   if(!ifs.is_open()) {
   std::cout << "codegen::checkConfigFile_ ERROR, unable to open "
   << config_file_name_ << std::endl;
   return;
   }
   }
   else {
   std::string alt_config_file_name = 
   source_header_name_.substr(0, source_header_name_.find_last_of("."));
   alt_config_file_name = alt_config_file_name + ".dg2";
   ifs.open(alt_config_file_name.c_str());
   if(!ifs.is_open()) {
   return;
   }
   }
   
  // reading the config file
  char buff[1024];
  while( ifs.good() ) {
  ifs.getline(buff, 1024);
  templateInstantiation ti;
  bool ret = ti.readFromConfigFile(buff);
  if(ret && (ti.name()==className)) tvec.push_back(ti);
  }
  
  ifs.close();
  }
*/

//bool codegen::parseCommandLine_(std::string const& buff, std::string const& rqst_className, TemplateInst_& ti)
//{
//  stringstream sstrm(buff.c_str());
//  std::string command, className;
//  std::vector<std::string> templateArgs_;
//  sstrm >> command >> className;
//  if( !sstrm.good() ||
//      (command!="register") || 
//      (className!=rqst_className) ) return false;
//  
//  // parse the template arguments
//  while(sstrm.good()) {
//    std::string nm,symt,realt;
//    std::string::size_type pos;
//    sstrm >> nm;
//    pos=nm.find_first_of("=");
//    symt = nm.substr(0,pos);
//    realt = nm.substr(pos+1,std::string::npos);
//    
//    if(symt=="" || realt=="") 
//      std::cout << "codegen::parseCommandLine_ ERROR parsing args"
//		<< std::endl;
//    return false;
//    ti.symbolicTypes_.push_back(symt);
//    ti.realTypes_.push_back(realt);
//  }
//
//  int i;
//  ti.fullName_ = className + "<";
//  if(ti.realTypes_.size()>0) 
//    ti.fullName_ = ti.fullName_ + ti.realTypes_[0];
//  for(i=0;i<(int)ti.realTypes_.size();i++)
//    ti.fullName_ = ti.fullName_ + ",";
//  ti.fullName_ = ti.fullName_ + ">";
//  
//  return true;
//}


void codegen::classPersisterDecl_(dict const& dict_)
{
  if(!dict_.isPersistent()) { return; }
  
  unsigned int i,j;
  std::string nm;
  
  // the includes 
  std::map<std::string,bool> types_;
  for(i=0;i<dict_.size();i++) {
    std::map<std::string,bool>::iterator it;
    it = types_.find(dict_.member(i).type);
    if(it == types_.end() ) {
      class_ofs_h_ << "#ifdef " << dict_.member(i).type << "__is__persistent" << std::endl
		   << "#include \"" << dict_.member(i).type << "__persister.h\"" << std::endl
		   << "#endif" << std::endl;
      types_[dict_.member(i).type]=true;
    }
  }
  class_ofs_h_ << std::endl << std::endl;
  
  // the read and write functions
  class_ofs_h_ << "template<class IOS" << dict_.templateSymbolicArgListDecl() << " >" << std::endl
	       << "void write(obj_output<IOS>& oo,"
	       << dict_.fullSymbolicName() << " const& p, const char* name=0)" << std::endl
	       << "{" << std::endl
	       << "  persister<" << dict_.fullSymbolicName() << ", IOS> pp(p);" << std::endl
	       << "  oo.write(pp);" << std::endl
	       << "}" << std::endl << std::endl;
  
  class_ofs_h_ << "template<class IOS" << dict_.templateSymbolicArgListDecl() << " >" << std::endl
	       << "void read(obj_input<IOS> const& oi,"
	       << dict_.fullSymbolicName() << "& p)" << std::endl
	       << "{" << std::endl
	       << "  persister<" << dict_.fullSymbolicName() << ",IOS> pp(p);" << std::endl
	       << "  oi.read(pp);" << std::endl
	       << "}" << std::endl << std::endl;

  // the << and >> operators
  class_ofs_h_ << "template<class IOS" << dict_.templateSymbolicArgListDecl() << " >" << std::endl
	       << "obj_output<IOS>& operator<<(obj_output<IOS>& oo,"
	       << dict_.fullSymbolicName() << " const& p)" << std::endl
	       << "{" << std::endl
	       << "  persister<" << dict_.fullSymbolicName() << ", IOS> pp(p);" << std::endl
	       << "  oo.write(pp);" << std::endl
	       << "  return oo;" << std::endl
	       << "}" << std::endl << std::endl;
  
  class_ofs_h_ << "template<class IOS" << dict_.templateSymbolicArgListDecl() << " >" << std::endl
	       << "obj_input<IOS> const& operator>>(obj_input<IOS> const& oi,"
	       << dict_.fullSymbolicName() << "& p)" << std::endl
	       << "{" << std::endl
	       << "  persister<" << dict_.fullSymbolicName() << ",IOS> pp(p);" << std::endl
	       << "  oi.read(pp);" << std::endl
	       << "  return oi;" << std::endl
	       << "}" << std::endl << std::endl;
  
  class_ofs_h_ << "template<class IOS" << dict_.templateSymbolicArgListDecl() << " >" << std::endl
	       << "class persister<" << dict_.fullSymbolicName() << ",IOS> : "
	       << "public handle<" << dict_.fullSymbolicName() << " >, "
	       << "public persister_base<IOS>  {" << std::endl;
  
  class_ofs_h_ << "public:" << std::endl;
  
  // constructors & destructor
  class_ofs_h_ << "  persister<" 
	       << dict_.fullSymbolicName() << ",IOS>() : handle<" 
	       << dict_.fullSymbolicName() << " >() {}" << std::endl;
  
  class_ofs_h_ << "  persister<" 
	       << dict_.fullSymbolicName() << ",IOS>("
	       << dict_.fullSymbolicName() << "& obj) : handle<"
	       << dict_.fullSymbolicName() << " >(&obj) {}" << std::endl;
    
  class_ofs_h_ << "  persister<" << dict_.fullSymbolicName() << ",IOS>("
	       << dict_.fullSymbolicName() << " const & obj) : handle<"
	       << dict_.fullSymbolicName() << " >(const_cast<"
	       << dict_.fullSymbolicName() << "*>(&obj)) {}// I know, that's UGLY." << std::endl;
  
  class_ofs_h_ << "  ~persister<" << dict_.fullSymbolicName() 
	       << ",IOS>() {}" << std::endl << std::endl;
  
  // interface methods
  
  // virtual unsigned int version() const { return <version>; }
  class_ofs_h_ << "  virtual unsigned int version() const { return " 
	       << dict_.version() << "; }" << std::endl;
  
  // virtual std::string name() const { return (std::string)<name>; }
  class_ofs_h_ << "  virtual std::string  name()    const { return (std::string)className_; }"
	       << std::endl;
  
  // virtual persister_base<IOS>* clone() const { return new persister<...,IOS>; }
  class_ofs_h_ << "  virtual persister_base<IOS>* clone() const { return new persister<"
	       << dict_.fullSymbolicName() << ",IOS>; }" << std::endl;
  
  // virtual void const*  get_object_addr() const { return (void*)obj_; }
  class_ofs_h_ << "  virtual void const* get_object_addr() const { return (void*)obj_; }"
	       << std::endl;
  
  // unsigned int size() const { return <size>; }
  class_ofs_h_ << "  unsigned int size()    const { return " 
	       << dict_.size() << ";}" << std::endl;
  
  // std::string name(unsigned int i) const { return (std::string)memberNames_[i]; }
  class_ofs_h_ << "  std::string  name(unsigned int i) const { return (std::string)memberNames_[i]; }" 
	       << std::endl;
  
  // std::string type(unsigned int i) const { return (std::string)memberTypes_[i]; }
  class_ofs_h_ << "  std::string  type(unsigned int i) const { return (std::string)memberTypes_[i]; }" 
	       << std::endl << std::endl;
  
  // 
  // now, the private section of the persister
  // 
  class_ofs_h_ << "private:" << std::endl;
  
  
  // the read_members() and write_members() methods 
  class_ofs_h_ << "  virtual void write_members(obj_output<IOS>& oo) const {"
	       << std::endl;
  for(i=0;i<dict_.baseList().size();i++) {
    nm = dict_.baseList()[i];
    class_ofs_h_ << "    {persister<" << nm << ",IOS> p(*(" << nm << "*)obj_);"
		 << "write(oo,p);}" << std::endl;
  }
  for(i=0;i<dict_.size();i++)
    if( dict_.member(i).isPersistent && dict_.member(i).arraySize.size()==0 )
      class_ofs_h_ << "    write(oo, obj_->" << dict_.member(i).name
		   << ", \"" << dict_.member(i).name << "\");" << std::endl;
    else if( dict_.member(i).isPersistent && dict_.member(i).arraySize.size()==1 )
      class_ofs_h_ << "    write(oo, obj_->" << dict_.member(i).name
		   << "," << dict_.member(i).arraySize[0] <<  ", \"" << dict_.member(i).name << "\");" << std::endl;
  
  class_ofs_h_ << "}" 
	       << std::endl << std::endl;
  
  class_ofs_h_ << "  virtual void read_members(obj_input<IOS> const& oi) {"
	       << std::endl;
  for(i=0;i<dict_.baseList().size();i++) {
    nm = dict_.baseList()[i];
    class_ofs_h_ << "     {persister<" << nm << ",IOS> p(*(" << nm << "*)obj_);"
		 << "read(oi,p);}" << std::endl;
  }
  for(i=0;i<dict_.size();i++) 
    if( dict_.member(i).isPersistent && dict_.member(i).arraySize.size()==0 )
      class_ofs_h_ << "    read(oi, obj_->" << dict_.member(i).name << ");" << endl;
    else if( dict_.member(i).isPersistent && dict_.member(i).arraySize.size()==1 )
      class_ofs_h_ << "    read(oi, obj_->" << dict_.member(i).name
		   << "," << dict_.member(i).arraySize[0] <<  ");" << std::endl;
  
  class_ofs_h_ << "}"
	       << std::endl << std::endl;
  
  // the friends
  class_ofs_h_ << "  template<class ZZ> friend class obj_input;"  << std::endl
	       << "  template<class ZZ> friend class obj_output;" 
	       << std::endl << std::endl;
  
  // the class description -- used by the schema evolution mechanism
  class_ofs_h_ << "  static const char* className_;" << std::endl
	       << "  static const char* memberNames_[];" << std::endl
	       << "  static const char* memberTypes_[];" << std::endl;
  
  // end of persister<OBJ,IOS>
  class_ofs_h_ << "};" << std::endl << std::endl;
  
  ofs_h_ << "#include \"" << dict_.name() << "__persister.h\"" << std::endl;
}


void codegen::classPersisterCode_(dict const& dict_)
{
  if(!dict_.isPersistent()) return;
  
  unsigned int i;
  
  ofs_cc_ << "const char* persister<" << dict_.fullRealName() << " ,xmlstream>::className_ = \""
	  << dict_.fullRealName() << "\";" << std::endl;
  
  ofs_cc_ << "const char* persister<" << dict_.fullRealName() << " ,xmlstream>::memberNames_[] = {";
  if(dict_.size()>0) {
    ofs_cc_ << "\""  << dict_.member(0).name  << "\"";
    for(i=1;i<dict_.size();i++)
      ofs_cc_ << ", \""  << dict_.member(i).name  << "\"";
  }
  ofs_cc_ << "};" << std::endl;
  
  ofs_cc_ << "const char* persister<" << dict_.fullRealName() << " ,xmlstream>::memberTypes_[] = {";
  if(dict_.size()>0) {
    unsigned int i;
    ofs_cc_ << "\""  << dict_.member(0).type  << "\"";
    for(i=1;i<dict_.size();i++)
      ofs_cc_ << ", \""  << dict_.member(i).type  << "\"";
  }
  ofs_cc_ << "};" << std::endl;
  
  ofs_cc_ << "type_registrar<" << dict_.fullRealName() << ",xmlstream> __" << dict_.name() << counter_++ << "__registration__;"
	  << std::endl;
  
  ofs_cc_ << std::endl;
}














#ifdef GARBAGE
  //  << "void write_members(obj_output<IOS>& oo) const {" << std::endl;
  
  
  unsigned int i;
  std::string nm;
  // first, we write the base classes
  for(i=0;i<(int)dict_.baseList().size();i++) {
    nm = dict_.baseList()[i];
    class_ofs_h_ << "{persister<" << nm << "> p(*("
		 << nm << "*)obj_);" 
		 << "oo.write(p);}" << std::endl;
  }
  for(i=0;i<dict_.size();i++) 
    if(dict_.member(i).isPersistent)
      class_ofs_h_ << "oo.write(obj_->"  << dict_.member(i).name  
		   << ", \""  << dict_.member(i).name  
		   << "\"" << ");" << std::endl;
  
  class_ofs_h_ << "}" << std::endl
	       << std::endl
	       << "template<class IOS>"
	       << "void read_members(obj_input<IOS> const& oi) {" << std::endl;
  
  for(i=0;i<(int)dict_.baseList().size();i++) {
    nm = dict_.baseList()[i];
    class_ofs_h_ << "{persister<" << nm << "> p(*("
		 << nm << "*)obj_);" 
		 << "oi.read(p);}" << std::endl;
  }
  for(i=0;i<dict_.size();i++)
    if(dict_.member(i).isPersistent)
      class_ofs_h_ << "oi.read(obj_->"  << dict_.member(i).name  
		   << ");" << std::endl;
  
  class_ofs_h_ << "}" << std::endl
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
	       << "persister<" << dict_.name() << " > pp(p);" << std::endl
	       << "oi.read(pp);" << std::endl
	       << "return oi;" << std::endl
	       << "};" << std::endl;
#endif
