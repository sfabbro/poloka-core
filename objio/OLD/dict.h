// -*- C++ -*-
// 
// file dict.h
// 
// 
#ifndef DICT_H
#define DICT_H

#include <string>
#include <vector>
#include <list>

#include <libxml/xmlreader.h>
#include <libxml/xmlwriter.h>


class dict {
public:
  
  struct persistentMember {
    persistentMember() 
      : name(""), type(0), typeName(""),
	scope(0), kind(""), isClean(true) {}
    ~persistentMember() {}
    
    std::string   name;
    unsigned int  type;
    std::string   typeName;
    unsigned int  scope;
    std::string   kind; // kind of scope...
    bool          isClean; // no raw pointers
    
    bool operator=(persistentMember const& pm) const {
      return (name==pm.name && typeName==pm.typeName);
    }
    
    void print() const {
      std::cout << typeName << " " << name << ";"
		<< std::endl;
    }
  };


  dict();
  ~dict();

  unsigned int                     size() const { return members_.size(); }
  
  unsigned int                     version() const { return version_; }
  unsigned int&                    version()       { return version_; }
  
  std::string const&               name() const { return name_; }
  std::string&                     name()       { return name_; }
  
  std::string const&               kind() const { return kind_; }
  std::string&                     kind()       { return kind_; }
  
  bool                             isTemplate() const { return isTemplate_; }
  bool&                            isTemplate()       { return isTemplate_; }
  
  bool                             isPersistent() const { return isPersistent_; }
  bool                             isPersistent()       { return isPersistent_; }
  
  std::list<std::string> const&    inheritanceList() const { return inheritanceList_; }
  std::list<std::string>&          inheritanceList()       { return inheritanceList_; }
  
  persistentMember const&          member(unsigned int i) const { return members_[i]; }
  void                             addPersistentMember(persistentMember const& pm) { members_.push_back(pm); }
  
  int                              writeXMLSchema(const std::string& filename) const;
  
  int                              readXMLSchema(const std::string& filename);
  
  int                              readFromSwigXMLStream(xmlTextReaderPtr);
  
  bool                             operator=(dict const&) const;
  
  void                             print(int verbosity) const;
  
  void                             clear();
  
  void                             test(xmlTextReaderPtr reader);
  
  
private:
  unsigned int    version_;
  std::string     name_;
  std::string     kind_;
  bool            isTemplate_;
  bool            isPersistent_;
  std::list<std::string> inheritanceList_;
  std::vector<persistentMember> members_;
  
  int  nextOpeningTag_(xmlTextReaderPtr reader, xmlChar const* name=0);
  int  nextClosingTag_(xmlTextReaderPtr, xmlChar const* name=0);
  int  nextElement_(xmlTextReaderPtr,xmlChar const* name=0);
  
  int  readClassAttributeList_(xmlTextReaderPtr);
  int  readClassBody_(xmlTextReaderPtr);
  int  readCDeclAttributeList_(xmlTextReaderPtr);
  int  readInheritanceList_(xmlTextReaderPtr);
  void decodeAttribute_(xmlTextReaderPtr, std::string& name, std::string& value);
  
  void printNodeInfo_(xmlTextReaderPtr reader); // debug function
  
  void        pad_(xmlTextWriterPtr writer, unsigned int depth) const;
  std::string tr_type_name_(std::string const&);
};


#endif

