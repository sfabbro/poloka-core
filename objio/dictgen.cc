/**
 *
 * dictgen.c 
 *
 * parse the swig-1.3 XML output. Identifies the class declarations
 * and attempts to generate persister and xml-schema descriptions
 * of the data structures.
 * 
 */

#include <unistd.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <list>
#include <vector>

#include <libxml/xmlreader.h>

using namespace std;



class PersistentClassDescription {
public:
  PersistentClassDescription() 
    : name_(""), kind_(""), isTemplate_(false), sz_(0) { }
  ~PersistentClassDescription() { }
  
  string const&          kind() const { return kind_; }
  string const&          name() const { return name_; }
  bool                   isTemplate() const { return isTemplate_; }
  vector<string> const&  inheritanceList() const { return baseList_; }
  
  unsigned int           size() const { return sz_; } 
  string const&          memberName(unsigned int i) const { return memberNames_[i]; };
  string const&          memberType(unsigned int i) const { return memberTypes_[i]; };
  
  int                    readFromSwigXMLStream(xmlTextReaderPtr);
  //  void           writeXMLSchema();
  
  void                   print() const;
  void                   test(xmlTextReaderPtr reader) const;
  
private:
  string name_;
  string kind_;
  bool   isTemplate_;
  unsigned int sz_;
  vector<string>  baseList_;
  vector<string>  memberNames_;
  vector<string>  memberTypes_;
  
  int    nextElement_(xmlTextReaderPtr, xmlChar const* name=0) const;
  
  int    readClassAttributeList_(xmlTextReaderPtr);
  int    readInheritanceList_(xmlTextReaderPtr);
  void   decodeAttribute_(xmlTextReaderPtr, string& name, string& value);
  int    readClassBody_(xmlTextReaderPtr);
  int    readCDeclAttributeList_(xmlTextReaderPtr);
  
  void   clear_() {
    name_="";
    kind_="";
    isTemplate_=false;
    sz_=0;
    baseList_.clear();
    memberNames_.clear();
    memberTypes_.clear();
  }
};


void PersistentClassDescription::test(xmlTextReaderPtr reader) const
{
  int ret;
  xmlChar* xmlc;
  do {
    ret = nextElement_(reader);
    xmlc = xmlTextReaderName(reader);
    cout << " ret=" << ret << ": " << (char*)xmlc << endl;
  } while(ret > 0 );
}


int PersistentClassDescription::readFromSwigXMLStream(xmlTextReaderPtr reader)
{
  clear_();
  if(!reader) return -1;
  int ret;
  xmlChar* name;
  
  ret = nextElement_(reader, (xmlChar*)"class");
  if(ret<=0) return ret;
  
  ret = readClassAttributeList_(reader);
  if(ret<=0) {
    std::cout << "[readFromSwigXMLStream] ERROR parsing the class attribute list." 
	      << std::endl;
    return ret;
  }
  
  ret = readClassBody_(reader);
  //  if(ret<=0) return ret;
  
  // closing tag...
  //  ret = nextElement_(reader, (xmlChar*)"class");
  //  if(ret!=2) 
  //    std::cout << "Ooops!" << std::endl;
  
  return 1;
}


int PersistentClassDescription::nextElement_(xmlTextReaderPtr reader, xmlChar const* name) const
{
  if(!reader) return -1;
  int ret, type, final_ret=-1;
  
  while(1) {
    ret = xmlTextReaderRead(reader);
    if(ret<=0) return ret;
    type = xmlTextReaderNodeType(reader);
    if(type == XML_READER_TYPE_ELEMENT ) final_ret=1;
    else if(type == XML_READER_TYPE_END_ELEMENT) final_ret=2;
    else continue;
    if(!name) return final_ret;
    xmlChar const* curName = xmlTextReaderConstName(reader);
    if( xmlStrcmp(curName, name) == 0 ) return final_ret;
  };
}


int PersistentClassDescription::readClassAttributeList_(xmlTextReaderPtr reader)
{
  int ret;
  std::string name_str, value_str;
  
  ret = nextElement_(reader);
  if(ret<=0) return ret;
  
  xmlChar const* eltname = xmlTextReaderConstName(reader);
  if( xmlStrcmp(eltname, (xmlChar*)"attributelist") != 0) {
    std::cout << "[readClassAttributeList_] ERROR tag <attributelist> not found"
	      << "(element name=" << (char*)eltname << std::endl;
    return -1;
  }
  
  while(1) {
    ret = nextElement_(reader);
    if(ret<=0) {
      std::cout << "[readClassAttributeList_] ERROR parsing the class attribute list"
		<< std::endl;
      return -1;
    }
    if(ret==2) 
      break;
    
    const xmlChar* name = xmlTextReaderConstName(reader);
    if( xmlStrcmp(name, (xmlChar*)"attribute") == 0 ) {
      decodeAttribute_(reader, name_str, value_str);
      if(name_str=="name") name_ = value_str;
      if(name_str=="kind") kind_ = value_str;
    }
    else if( xmlStrcmp(name, (xmlChar*)"baselist") == 0 )
      readInheritanceList_(reader);
  }
  
  return 1;
}


int PersistentClassDescription::readCDeclAttributeList_(xmlTextReaderPtr reader)
{
  int ret;
  std::string name_str, value_str;
  std::string el_name="", el_type="";
  bool skip = false;
  
  ret = nextElement_(reader);
  if(ret<=0) return ret;
  
  xmlChar const* eltname = xmlTextReaderConstName(reader);
  if( xmlStrcmp(eltname, (xmlChar*)"attributelist") != 0) {
    std::cout << "[readCDeclAttributeList_] ERROR tag <attributelist> not found"
	      << "(element name=" << (char*)eltname << std::endl;
    return -1;
  }
  
  // FIXME: throw away thse depth tests and 
  // implement a real nextSibling...
  int cur_depth, depth=xmlTextReaderDepth(reader);
  
  while(1) {
    ret = nextElement_(reader);
    if(ret<=0) {
      std::cout << "[readCDeclAttributeList_] ERROR parsing CDecl attribute list"
		<< std::endl; 
      return -1;
    }
    cur_depth = xmlTextReaderDepth(reader);
    if(cur_depth>(depth+1)) continue;
    
    if(ret==2)
      break;
    
    const xmlChar* name = xmlTextReaderConstName(reader);
    if( xmlStrcmp(name, (xmlChar*)"attribute") == 0 ) {
      decodeAttribute_(reader, name_str, value_str);
      if(name_str=="name") el_name = value_str;
      if(name_str=="type") el_type = value_str;
      if(name_str=="code") skip = true;
    }
  }
  
  if(el_name=="" || el_type=="" || skip) return 1;
  
  memberNames_.push_back(el_name);
  memberTypes_.push_back(el_type);
  sz_++;
  
  return 1;
}


int PersistentClassDescription::readInheritanceList_(xmlTextReaderPtr reader)
{
  int ret;
  
  while(1) {
    ret = nextElement_(reader);
    if(ret<=0) {
      std::cout << "[readInheritanceList_] ERROR while parsing class Inheritance list"
		<< std::endl;
      return -1;
    }
    if(ret==2) break;
    
    xmlChar* name = xmlTextReaderGetAttribute(reader, (xmlChar*)"name");
    baseList_.push_back((char*)name);
    xmlFree(name);
  }
  
  return 1;
}


int PersistentClassDescription::readClassBody_(xmlTextReaderPtr reader)
{
  int ret;
  xmlChar* elementName;
  
  while(1) {
    ret = nextElement_(reader);
    if(ret<=0) return ret;
    elementName = xmlTextReaderName(reader);
    if(ret==2 && xmlStrcmp(elementName, (xmlChar*)"class")==0 )
      return 1;
    if( xmlStrcmp(elementName, (xmlChar*)"cdecl")!=0 ) 
      continue;
    if(ret!=1) continue;
    readCDeclAttributeList_(reader);
  }
}



void PersistentClassDescription::decodeAttribute_(xmlTextReaderPtr reader, string& name, string& value)
{
  xmlChar* name_tmp = xmlTextReaderGetAttribute(reader, (xmlChar*)"name");
  xmlChar* value_tmp = xmlTextReaderGetAttribute(reader, (xmlChar*)"value");
  
  name = (char*)name_tmp;
  value = (char*)value_tmp;
  
  xmlFree(name_tmp);
  xmlFree(value_tmp);
}


void PersistentClassDescription::print() const
{
  unsigned int i;
  
  if(isTemplate_) cout << "template<class T>";
  cout << kind_ << " " << name_;
  if(baseList_.size()>0) {
    cout << " : public ";
    for(i=0;i<baseList_.size();i++) {
      cout << baseList_[i];
      if(i<(baseList_.size()-1)) cout << ", ";
    }
  }
  cout << ";" << endl ;
  
  for(i=0;i<sz_;i++)
    cout << " [" << setw(8) << memberTypes_[i] << "  " << memberNames_[i] << "]" << endl;
  
  
}






void usage()
{
  fprintf(stderr, "usage: dictgen [OPTIONS] <filename>\n");
  fprintf(stderr, " where OPTIONS can be:\n");
  fprintf(stderr, "     -v                     increase verbosity level\n");
  fprintf(stderr, "     -o <outfile_base_name> set the output base name\n");
  fprintf(stderr, "     <filename>             input file name (stdin if -)\n");
  exit(-1);
}



int main(int argc, char** argv)
{
  char c;
  int verb=0, ret;
  std::string filename="";
  std::string outfileBaseName="";
  xmlTextReaderPtr reader;
  
  while( (c=getopt(argc, argv, "vo:")) != -1 )
    switch(c) {
    case 'v':
      verb++;
      break;
    case 'o':
      outfileBaseName=optarg;
      break;
    default:
      usage();
    }
  if(optind>=argc) usage();
  filename=argv[optind];
  if(outfileBaseName=="") 
    outfileBaseName="dict_"+filename;
  if(verb)
    std::cout << "dictgen " << filename << " --> "
	      << outfileBaseName << ".[h|cc]"<< std::endl;
  
  // open the reader ... 
  reader = xmlNewTextReaderFilename(filename.c_str());
  if( !reader ) {
    std::cerr << "dictgen: unable to open file. Aborting."
	      << std::endl;
    exit(-1);
  }
  
  PersistentClassDescription pdd;
  
  pdd.readFromSwigXMLStream(reader);
  pdd.print();
  pdd.readFromSwigXMLStream(reader);
  pdd.print();
  pdd.readFromSwigXMLStream(reader);
  pdd.print();
  
  xmlFreeTextReader(reader);
}

