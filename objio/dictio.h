// -*- C++ -*-
// 
// file dictio.h
// 
// 
#ifndef DICTIO_H
#define DICTIO_H

#include <list>
#include <vector>
#include <map>

#include "dict.h"

template<class T> class persister;


template<class IOS>
class dict_output {
public:
  dict_output(const std::string& name, int compression=0) : stream_(name, compression) {}
  dict_output(IOS& stream) : stream_(stream) {}
  ~dict_output() { close(); }
  
  void         open(std::string const& filename, int compression) const {
    // FIXME: close the stream first
    stream_.open(filename, compression);
  }
  
  void         close() { stream_.close(); }
  
  template<class T>
  void    write(persister<T>& pt) {
    stream_.start_dict(pt.name(), pt.version());
    //    pt.write_members(*this);
    stream_.end_dict();
  }
  
  virtual void   write(int1 v,   const char* name) { stream_.write_member(name, "int1"); }
  virtual void   write(uint1 v,  const char* name) { stream_.write_member(name, "uint1"); }
  virtual void   write(int2 v,   const char* name) { stream_.write_member(name, "int2"); }
  virtual void   write(uint2 v,  const char* name) { stream_.write_member(name, "uint2"); }
  virtual void   write(int4 v,   const char* name) { stream_.write_member(name, "int4"); }
  virtual void   write(uint4 v,  const char* name) { stream_.write_member(name, "uint4"); }
  virtual void   write(int8 v,   const char* name) { stream_.write_member(name, "int8"); }
  virtual void   write(uint8 v,  const char* name) { stream_.write_member(name, "uint8"); }
  virtual void   write(float4 v, const char* name) { stream_.write_member(name, "float4"); }
  virtual void   write(float8 v, const char* name) { stream_.write_member(name, "float8"); }  
  virtual void   write(std::string const& v, const char* name) { stream_.write_member(name, "string"); }
  
  template<class T>
  void           write(std::list<T> const& l, const char* name) {
    persister<T> pp;
    std::stringstream sstrm;
    sstrm << "collection<" << pp.type << "[" << pp.version << "]>";
    stream_.write_member(name, sstrm.str().c_str());
  }
  
  template<class T>
  void           write(std::vector<T> const& v, const char* name) {
    persister<T> pp;
    std::stringstream sstrm;
    sstrm << "collection<" << pp.type << "[" << pp.version << "]>";
    stream_.write_member(name, sstrm.str().c_str());
  }
  
  template<class T, class U>
  void           write(std::map<T,U> const& m, const char* name=0) {
    persister<T> pp;
    std::stringstream sstrm;
    sstrm << "collection<" << pp.type << "[" << pp.version << "]>";
    stream_.write_member(name, sstrm.str().c_str());
  }
  
  
  template<class T, class U>
  void           write(std::pair<T,U> const& p, const char* name=0) {
    persister<T> pt;
    persister<U> pu;
    std::stringstream sstrm;
    sstrm << "collection.pair<" 
	  << pt.type() << "[" << pt.version() << "], "
	  << pu.type() << "[" << pu.version() << "]>";
    stream_.write_member(name, sstrm.str().c_str());
  }
  
  
private:
  IOS stream_;
};


template<class IOS>
class dict_input {
public:
  dict_input(const std::string& name) : stream_(name) {}
  ~dict_input() { }
  
  void       open(std::string const& filename) const {
    stream_.open(filename);
  }
  
  void         close() { stream_.close(); }
  //  inline void  read(dict&) const;
  
private:
  IOS stream_;
};



/////////////////////////////////////// INLINED methods //////////////////////////////////
/* 
   template<class IOS>
   void dict_output<IOS>::write(dict const& d)
   {
   if(!d.isPersistent()) return;
   stream_.start_dict(d.name(), d.kind(), d.size(), d.version());
   stream_.start_collection(d.inheritanceList().size(), "inheritanceList");
   std::list<std::string>::const_iterator ilit;
   for(ilit=d.inheritanceList().begin();ilit!=d.inheritanceList().end();ilit++)
   stream_.write(*ilit);
   stream_.end_collection();
   stream_.start_collection(d.size(), "memberList");
   unsigned int i;
   for(i=0;i<d.size();i++) {
   stream_.write(d.member(i).name, "name");
   stream_.write(d.member(i).type, "type");
   stream_.write(d.member(i).typeName, "typeName");
   stream_.write(d.member(i).scope, "scope");
   stream_.write(d.member(i).kind, "kind");
   stream_.write(d.member(i).isClean, "isClean");
   }
   stream_.end_collection();
   stream_.end_dict();
   }
   
   
   
   template<class IOS>
   void dict_input<IOS>::read(dict& d) const
   {
   d.clear();
   std::string name, kind, baseClass;
   unsigned int i, size, version;
   
   stream_.start_dict(name, kind, size, version);
   d.name()=name; d.kind()=kind; d.version()=version;
   
   stream_.start_collection(size);
   std::list<std::string>::const_iterator ilit;
   for(i=0;i<size;i++) {
   stream_.read(baseClass);
   d.inheritanceList().push_back(baseClass);
   }
   stream_.end_collection();
   
   stream_.start_collection(size);
   for(i=0;i<size;i++) {
   //    dict::persistentMember pm;
   //    stream_.read(pm.name);
   //    stream_.read(pm.type);
   //    stream_.read(pm.typeName);
   //    stream_.read(pm.scope);
   //    stream_.read(pm.kind);
   //    stream_.read(pm.isClean);
   //    d.addPersistentMember(pm);
   }
   
   stream_.end_collection();
   stream_.end_dict();
   }
*/

#endif
