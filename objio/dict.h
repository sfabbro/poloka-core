// -*- C++ -*-
// 
// file dict.h
// 
// 
#ifndef DICT_H
#define DICT_H

#include <string>
#include <vector>

class dict_writer_base;
class dict_reader_base;


class dict {
public:

  struct member_t {
    member_t();
    ~member_t() { clear(); }
    
    std::string name;
    std::string type;
    std::vector<unsigned int> arraySize;
    bool        isPointer;
    bool        isStatic;
    bool        isTemplate;
    bool        complexType;
    bool        isPersistent;
    bool        isReference;
    
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
  bool                             isTemplate()   const { return isTemplate_; }
  bool                             isPersistent() const { return isPersistent_; }
  std::string const&               kind()         const { return kind_; }
  std::vector<std::string> const&  baseList() const { return baseList_; }
  
  std::string const&               memberName(unsigned int i) const { return memberList_[i].name; }
  std::string const&               memberType(unsigned int i) const { return memberList_[i].type; }
  dict::member_t                   member(unsigned int i) const { return memberList_[i]; }
  
  bool                             operator==(dict const&) const;
  
  void                             write(dict_writer_base&) const;
  void                             read(dict_reader_base const&);
  
  void                             print(int verbosity) const;
  

protected:
  unsigned int version_;
  std::string  name_;
  std::string kind_;
  bool isTemplate_;
  bool isPersistent_;
  std::vector<std::string> baseList_;
  std::vector<member_t>    memberList_;
  
  void         clear_();
  
  friend class dict_reader_base;
  friend class dict_writer_base;

};


#endif

