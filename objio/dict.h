// -*- C++ -*-
// 
// file dict.h
// 
// The dict class is used to describe the structure
// of the persistent parts of a C++ class. A dict can 
// be read by the dict_reader classes such as swig_dict_reader.
// 
// NOTE: this class is no more than a quick and dirty trick and should
// be rewritten cleanly in a near future. We should implement a
// CppType class, to manage cleanly the C++ type specifications and
// the template instantiation mechanism. A PersistentMember would then
// be a {CppType,CppName} struct, and a dict would just be a class like:
// 
// Dict {
//  list<string>     templateArgs;
//  string           name();
//  unsigned int     version();
//  PersistentMember m[];
// };
// 
#ifndef DICT_H
#define DICT_H

#include <string>
#include <vector>

class dict_writer_base;
class dict_reader_base;

class templateInstantiation;

class dict {
public:

  struct member_t {
    member_t();
    ~member_t() { clear(); }
    
    std::string name;           // the member name (ex. x, y, flux)
    std::string type;           // the full type name (ex. double*, float**, Star[], std::vector<T>)
    std::string baseType;       // the base type (double,float,Star)
    
    std::vector<unsigned int> arraySize; // useless (arrays are not persistent. use STL types instead)
    
    bool        isPointer;      // true if pointer type (raw pointer members are not persistent!)
    bool        isStatic;       // true if static member (static members are not persistent!)
    bool        isTemplate;     // useless ?
    bool        complexType;    // useless ?
    bool        isReference;    // true if member is a reference. (references are not persistent!)
    bool        isPersistent;   // true if the member is persistent
    
    void        print() const;
    void        update();
    member_t&   operator=(member_t const& m) { copy(m); return *this; }
    void        copy(member_t const&);
    void        clear();
  };
  
public:
  dict();
  ~dict();
  
  unsigned int                     size()         const { return memberList_.size(); }
  unsigned int                     version()      const { return version_; }
  std::string const&               name()         const { return name_; }
  std::string const&               fullSymbolicName() const { return symName_; }
  std::string const&               fullRealName() const { return realName_; }
  bool                             isTemplate()   const { return isTemplate_; }
  bool                             isPersistent() const { return isPersistent_; }
  
  std::string const&               kind()         const { return kind_; }
  std::vector<std::string> const&  baseList()     const { return baseList_; }
  
  std::string const&               memberName(unsigned int i) const { return memberList_[i].name; }
  std::string const&               memberType(unsigned int i) const { return memberList_[i].type; }
  dict::member_t                   member(unsigned int i) const { return memberList_[i]; }
  int                              nTemplateArgs() const { return (int)symbolicTypes_.size(); }
  std::string                      templateSymbolicArg(int i) const { return symbolicTypes_[i]; }
  std::string                      templateRealArg(int i) const { return realTypes_[i]; }
  
  bool                             operator==(dict const&) const;
  
  void                             update(templateInstantiation const&);
  void                             instantiate(templateInstantiation const&);
  
  void                             write(dict_writer_base&) const;
  void                             read(dict_reader_base const&);
  
  void                             print(int verbosity) const;
  
  // FIXME: nrl -- 03/2004 put this in the code generator
  std::string                      templateSymbolicArgList() const;
  std::string                      templateSymbolicArgListDecl() const;
  
protected:
  unsigned int version_;
  std::string  name_;
  std::string  realName_,symName_;
  std::string kind_;
  bool isTemplate_;
  bool isPersistent_;
  std::vector<std::string> baseList_;
  std::vector<member_t>    memberList_;
  std::vector<std::string> symbolicTypes_;
  std::vector<std::string> realTypes_;
  
  void         clear_();
  
  friend class dict_reader_base;
  friend class dict_writer_base;
  
};


class templateInstantiation {
public:
  templateInstantiation();
  ~templateInstantiation();
  
  std::string    name() const { return templateClassName_; }
  int            nTemplateArgs() const { return (int)symbolicTypes_.size(); }
  

  std::string    fullSymbolicName() const;
  std::string    fullRealName() const;
  bool           hasSymbolicType(std::string const&) const;
  std::string    realType(std::string const& symtype) const;
  std::string    instantiateType(std::string const& sym) const;
  
  std::string    templateSymbolicArg(int i) const { return symbolicTypes_[i]; }
  std::string    templateRealArg(int i) const { return realTypes_[i]; }
  
  bool           readFromConfigFile(std::string const&);
  
  void           print() const;
  
  

  std::string templateClassName_;
  std::vector<std::string> symbolicTypes_;
  std::vector<std::string> realTypes_;

private:

  void clear_();
  void splitType_(std::string const&,std::vector<std::string>&) const;
  std::string rebuildType_(std::vector<std::string> const&) const;
};


#endif

