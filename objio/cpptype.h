// -*- C++ -*-
// 
// \file cpptype.h Simple class to hold the type of
// a class member.
// 
// 
// POTENTIAL PROBLEMS:
//   * the CppType does not know whether it is fully instantiated or not...
#ifndef CPPTYPE_H
#define CPPTYPE_H

#include <string>
#include <vector>


class CppTemplateInstance;


class CppType {
public:
  CppType();
  CppType(std::string const& decl) { readFromCppTypeDecl(decl); }
  ~CppType();
  
  //! the type name. For example: Point, Star, StarList<T>. 
  std::string        cppTypeName() const { return cppTypeName_; }
  
  //! the symbolic name (if the class is a template)
  std::string        symbolicCppTypeName() const { return symbolicCppTypeName_; }
  
  //! the base type (without the type modifiers, like *,&,const,[] etc...)
  std::string        baseCppTypeName() const { return baseCppTypeName_; }
  
  //! the symbolic class name, without the type modifiers
  std::string        symbolicBaseCppTypeName() const { return symbolicBaseCppTypeName_; }
  
  //! return the number of template args. return 0 if the class is not template
  int                nTemplateArgs() const { return templateArgs_.size(); }
  
  //! return the class template argument number i
  std::string        templateArg(int i) const { return templateArgs_[i]; }
  
  //! return the class template arg vec
  std::vector<std::string>  templateArgVec() const { return templateArgs_; }
  
  //! true if the class is a template classa
  bool               isTemplate() const { return isTemplate_; /* return templateArgs_.size()!=0; */ }
  
  //! true if the template is instantiated (StarList<SEStar>)
  bool               wasInstantiated() const { return wasInstantiated_; }
  
  //! true if the type is not an atomic type (int,double,...)
  bool               isAComplexType() const { return isAComplexType_; }
  
  //! true if the type is a reference (--> not persistent)
  bool               isAReference() const { return isAReference_; }
  
  //  //! true if the the type is a STL container (persistent)
  //  bool               isASTLContainer() const { return isASTLContainer_; }
  
  //! true if there is a const modifier somewhere in the type name (-> not persistent)
  bool               isConst() const { return isConst_; }
  
  //! true if this is a static type (-> not persistent)
  bool               isStatic() const { return isStatic_; }
  
  //! true if this is a raw pointer (-> not persistent)
  bool               isAPointer() const { return isAPointer_; }
  
  //! true if the type is persistent
  bool               isPersistent() const;
  
  //! instantiate the template class
  CppType            instantiate(CppTemplateInstance const&) const;
  
  //! copy the object 
  void               copy(CppType const&);
  
  //! copy operator --just calls copy()
  CppType&           operator=(CppType const& t) { copy(t); return *this; }
  
  //! initialize the object from the C++ type declaration
  void               readFromCppTypeDecl(std::string const&);
  
  //! print out the structure
  void               print() const;
  
  //! split the string str into tokens. if keep_del is true,
  //! the delimiters are kept in the tokens output vector
  static void         tokenize(std::string const& str,
			       std::vector<std::string>& tokens,
			       std::string const& del, bool keep_del);
  
  static std::string  buildTypeString(std::vector<std::string> const&);
  static std::string  reformatTypeString(std::string const&);
  static std::string  buildCleanTypeString(std::vector<std::string> const&);
  
  void                clear();
  
  
protected:
  std::string cppTypeName_;
  std::string symbolicCppTypeName_;
  std::string baseCppTypeName_;
  std::string symbolicBaseCppTypeName_;
  
  bool isConst_;
  bool isStatic_;
  bool isAPointer_;
  //  bool isASTLContainer_;
  bool isAComplexType_; // remove this?
  bool isAReference_;
  bool isTemplate_;
  bool wasInstantiated_; // means, the instantiate() function was called...
  
  std::vector<std::string> templateArgs_;  // all the template args
  std::vector<std::string> tok_;
  std::vector<std::string> baseTok_;
  
  static void  cleanTokens_(std::vector<std::string>&);
};




class CppTemplateInstance {
public:
  CppTemplateInstance() {}
  ~CppTemplateInstance() {}
  
  int                size() const { return sym_.size(); }
  
  std::string        realName(int i) const { return real_[i]; }
  std::string        realName(std::string const& sym) const;
  bool               hasRealName(std::string const& real) const { return find_(real_,real)>=0; }
  
  std::string        symName(int i) const { return sym_[i]; }
  std::string        symName(std::string const& real) const;
  bool               hasSymName(std::string const& sym) const { return find_(sym_,sym)>=0; }
  
  void               addInstance(std::string const& sym, std::string const& real);
  //  void               readFromHeaderSpec(std::string const&);
  
  void               print() const;
  
private:
  std::vector<std::string> sym_;
  std::vector<std::string> real_;
  
  void  clear_();
  int   find_(std::vector<std::string> const& v, std::string const& name) const;
};


#endif

