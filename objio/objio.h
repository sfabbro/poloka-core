// -*- C++ -*-
// $Id: objio.h,v 1.16 2004/03/16 17:31:22 guy Exp $
// 
// \file objio.h
// 
// Last modified: $Date: 2004/03/16 17:31:22 $
// by:            $Author: guy $
// 
#ifndef OBJIO_H
#define OBJIO_H

#include <typeinfo>
#include <string>
#include <vector>
#include <list>
#include <map>

#include "toadtypes.h"
#include "objio_defs.h"

#include "countedref.h"

using namespace std;

template<class IOS> class typemgr;
template<class IOS> class persister_base;
template<class T, class IOS> class persister;


template<class IOS>
class obj_output { //: public obj_output_base {
public:
  obj_output(const std::string& name) : stream_(name,OBJIO_WRITE,0) {}
  ~obj_output() { }

  void         open(std::string const& filename, int compression=0) const {
    stream_.open(filename,OBJIO_WRITE,compression);
  }

  void         close() { stream_.close(); }
  
  template<class T>
  void   write(persister<T,IOS> const& p, const char* name=0) {
    stream_.write_start_object_tag(p.name(), p.version());
    p.write_members(*this);
    stream_.write_end_object_tag();
  }
  
  void   write(persister_base<IOS> const* p, const char* name=0) {
    stream_.write_start_object_tag(p->name(), p->version());
    p->write_members(*this);
    stream_.write_end_object_tag();
  }
  
  virtual void   write(int1 const& v,   const char* name=0) { stream_.write(v, name); }
  virtual void   write(uint1 const& v,  const char* name=0) { stream_.write(v, name); }
  virtual void   write(int2 const& v,   const char* name=0) { stream_.write(v, name); }
  virtual void   write(uint2 const& v,  const char* name=0) { stream_.write(v, name); }
  virtual void   write(int4 const& v,   const char* name=0) { stream_.write(v, name); }
  virtual void   write(uint4 const& v,  const char* name=0) { stream_.write(v, name); }
  virtual void   write(int8 const& v,   const char* name=0) { stream_.write(v, name); }
  virtual void   write(uint8 const& v,  const char* name=0) { stream_.write(v, name); }
  virtual void   write(float4 const& v, const char* name=0) { stream_.write(v, name); }
  virtual void   write(float8 const& v, const char* name=0) { stream_.write(v, name); }
  virtual void   write(const std::string& v, const char* name=0) { stream_.write(v, name); }
  
  
  // tabs     //    stream_.write(t,sz,*this);
  template<class T>
  void           write(T const* t, unsigned long sz, const char* name=0) {
    stream_.write_start_collection_tag(sz,name);
    int i;
    for(i=0;i<sz;i++)
      write(*(t+i));
    stream_.write_end_collection_tag();
  }
  
  // raw pointers 
  template<class T>
  void           write(T const* t, const char* name=0) {
    bool w = check_address_((void const*)t); // did we write this object already ?
    persister<T,IOS>* b = (persister<T,IOS>*)typemgr<IOS>::getPersister_rtti(typeid(*t).name());
    b->set_object(t);
    stream_.write_start_raw_pointer_tag(name);
    if(!w) {
      write((persister_base<IOS>*)b);
    }
    else 
      stream_.write(t);
    stream_.write_end_raw_pointer_tag();
  }
  
  template<class T>
  void           write(std::list<T> const& l, const char* name=0) {
    stream_.write_start_collection_tag(l.size(), name);
    typename std::list<T>::const_iterator it;
    for(it=l.begin();it!=l.end();it++) {
      persister<T,IOS> p(*it);
      write(p);
      // *this << *it;
    }
    stream_.write_end_collection_tag();
  }

  template<class T>
  void           write(std::vector<T> const& v, const char* name=0) {
    stream_.write_start_collection_tag(v.size(), name);
    typename std::vector<T>::const_iterator it;
    for(it=v.begin();it!=v.end();it++) {
      persister<T,IOS> p(*it);
      write(p);
      //*this << *it;
    }
    stream_.write_end_collection_tag();
  }

  template<class T, class U>
  void           write( std::map<T,U> const& m, const char* name=0) {
    stream_.write_start_collection_tag(m.size(), name);
    typename std::map<T,U>::const_iterator it;
    for(it=m.begin();it!=m.end();it++) {
      // this time, we must use write().
      // Because we know what we are writing out...
      
      const  std::pair<T,U>&  pairTU = *it;
      write(pairTU);
      //write(*it);
    }
    stream_.write_end_collection_tag();
  }
  
  
  template<class T, class U>
  void           write( std::pair<T,U> const& p, const char* name=0) {
    stream_.write_start_collection_tag(2, name);
    persister<T,IOS> fp(p.first);
    persister<U,IOS> sp(p.second);
    write(fp);
    write(sp);
    stream_.write_end_collection_tag();
  }
  
  template<class T>
  void         write(CountedRef<T> const& r, const char* name=0) {
    T const* pt = &(*r);
    bool w = check_address_((void const*)(pt)); // did we write this object already ?
    persister<T,IOS>* b = (persister<T,IOS>*)typemgr<IOS>::getPersister_rtti(typeid(*pt).name());
    b->set_object(pt);
    stream_.write_start_reference_tag(name, b->get_object_addr());
    if(!w) {
      write((persister_base<IOS>*)b);
    }
    else 
      stream_.write(pt);
    stream_.write_end_reference_tag();
  }
  
  void         write(RefCount const& cr, const char* name=0) {
    stream_.write_start_object_tag("RefCount",0);
    stream_.write_end_object_tag();
  }
  
  //  template<class T>
  //    void           write(handle<T> const& h) {
  //    stream_.write(h,*this);
  //  }

  //  virtual void   write_streamer(streamer_base const* s) {
  //    stream_.write_streamer(s);
  //  }
  
private:
  IOS stream_;
  mutable std::map<void const*,bool> addr_;
  
  bool check_address_(void const* t) const {
    if(!t) return false;
    std::map<void const*,bool>::const_iterator it;
    it = addr_.find(t);
    if(it==addr_.end()) {
      addr_[t]=true;
      return false;
    }
    return true;
  }
  
  template<class T, class U> friend class persister;
};



#define define_output_operator(type)                       \
template<class IOS>                                        \
obj_output<IOS>& operator<<(obj_output<IOS>& oo, type const& v)   \
{                                                          \
  oo.write(v); return oo;                                  \
}                                                          \

define_output_operator(int1)
define_output_operator(uint1)
define_output_operator(int2)
define_output_operator(uint2)
define_output_operator(int4)
define_output_operator(uint4)
define_output_operator(int8)
define_output_operator(uint8)
define_output_operator(float4)
define_output_operator(float8)
define_output_operator(std::string)


#define define_output_function(type)                       \
template<class IOS>                                        \
void write(obj_output<IOS>& oo, type const& v, const char* name=0) \
{                                                          \
  oo.write(v,name);                                        \
}                                                          \

define_output_function(int1)
define_output_function(uint1)
define_output_function(int2)
define_output_function(uint2)
define_output_function(int4)
define_output_function(uint4)
define_output_function(int8)
define_output_function(uint8)
define_output_function(float4)
define_output_function(float8)
define_output_function(std::string)



template<class IOS, class T>
obj_output<IOS>& operator<<(obj_output<IOS>& oo, std::list<T> const& l)
{
  oo.write(l); return oo;
}

template<class IOS, class T>
obj_output<IOS>& operator<<(obj_output<IOS>& oo, std::vector<T> const& v)
{
  oo.write(v); return oo;
}

template<class IOS, class T, class U>
obj_output<IOS>& operator<<(obj_output<IOS>& oo, std::map<T,U> const& m)
{
  oo.write(m); return oo;
}

template<class IOS, class T>
obj_output<IOS>& operator<<(obj_output<IOS>& oo, CountedRef<T> const& cr)
{
  oo.write(cr); return oo;
}

template<class IOS>
obj_output<IOS>& operator<<(obj_output<IOS>& oo, RefCount const& cr)
{
  oo.write(cr); return oo;
}


template<class IOS, class T>
void write(obj_output<IOS>& oo, T const* p, unsigned long sz, const char* name=0)
{
  oo.write(p,sz,name);
}

template<class IOS, class T>
void write(obj_output<IOS>& oo, persister<T,IOS> const& p, const char* name=0)
{
  oo.write(p,name);
}

template<class IOS, class T>
void write(obj_output<IOS>& oo, persister_base<IOS> const* p, const char* name=0)
{
  oo.write(p,name);
}

template<class IOS, class T>
void write(obj_output<IOS>& oo, std::list<T> const& l, const char* name=0)
{
  oo.write(l,name);
}

template<class IOS, class T>
void write(obj_output<IOS>& oo, std::vector<T> const& v, const char* name=0)
{
  oo.write(v,name);
}

template<class IOS, class T, class U>
void write(obj_output<IOS>& oo, std::map<T,U> const& m, const char* name=0)
{
  oo.write(m,name);
}

template<class IOS, class T>
void write(obj_output<IOS>& oo, CountedRef<T> const& cr, const char* name=0)
{
  oo.write(cr,name);
}

template<class IOS>
void write(obj_output<IOS>& oo, RefCount const& cr, const char* name=0)
{
  oo.write(cr,name);
}

///////////////////////// OBJ_INPUT /////////////////////////

template<class IOS>
class obj_input { //: public obj_input_base {
public:
  obj_input(const std::string& name) : stream_(name,OBJIO_READ,0) {}
  ~obj_input() {}
  
  void         open(std::string const& filename) {
    stream_.open(filename,OBJIO_READ,0);
  }
  
  void         close() { stream_.close(); }
  
  //  template<class T>
  //  virtual persister_base<IOS>*  read() const;

  //  virtual persister_base<IOS>* read() const { return 0; }

  template<class T>
  void         read(persister<T,IOS>& p) const {
    unsigned int version;
    std::string name;
    stream_.read_start_object_tag(name, version); // may throw an exc. here
    p.read_members(*this);
    stream_.read_end_object_tag();
  }

  virtual void read(int1& v)   const { stream_.read(v); }
  virtual void read(uint1& v)  const { stream_.read(v); }
  virtual void read(int2& v)   const { stream_.read(v); }
  virtual void read(uint2& v)  const { stream_.read(v); }
  virtual void read(int4& v)   const { stream_.read(v); }
  virtual void read(uint4& v)  const { stream_.read(v); }
  virtual void read(int8& v)   const { stream_.read(v); }
  virtual void read(uint8& v)  const { stream_.read(v); }
  virtual void read(float4& v) const { stream_.read(v); }
  virtual void read(float8& v) const { stream_.read(v); }
  virtual void read(std::string& v) const { stream_.read(v); }

  virtual void skip() const { stream_.skip(); }

  template<class T>
  void         read(T*& t, unsigned long& sz) const {
    int i;
    unsigned int sz_tmp;
    T val;
    stream_.read_start_collection_tag(sz_tmp);
    sz_tmp = sz; // FIXME
    for(i=0;i<sz;i++) {
      *this >> t[i];
    }
    stream_.read_end_collection_tag();
  }

  template<class T>
  void         read(std::list<T>& l) const {
    unsigned int i, sz;
    T val;
    stream_.read_start_collection_tag(sz);
    l.clear(); // we reset the list
    for(i=0;i<sz;i++) {
      persister<T,IOS> p(val);
      read(p);
      //*this >> val;
      l.push_back(val);
    }
    stream_.read_end_collection_tag();
  }

  template<class T>
  void         read(std::vector<T>& v) const {
    unsigned int i, sz;
    T val;
    stream_.read_start_collection_tag(sz);
    v.clear();
    for(i=0;i<sz;i++) {
      persister<T,IOS> p(val);
      read(p);
      v.push_back(val);
    }
    stream_.read_end_collection_tag();
  }
  
  template<class T, class U>
  void         read(std::map<T,U>& m) const {
    unsigned int i, sz;
    std::pair<T,U> val;
    stream_.read_start_collection_tag(sz);
    m.clear();
    for(i=0;i<sz;i++) {
      read(val);
      m[val.first]=val.second;
    }
    stream_.read_end_collection_tag();
  }
  
  template<class T, class U>
  void         read(std::pair<T,U>& p) const {
    unsigned int sz;
    T v1; U v2;
    stream_.read_start_collection_tag(sz);
    persister<T,IOS> p1(v1);
    persister<U,IOS> p2(v2);
    read(p1);
    read(p2);
    //*this >> v1; 
    //*this >> v2;
    p.first=v1;
    p.second=v2;
    stream_.read_end_collection_tag();
  }
  
 template<class T>
  T* read_next_object() const{
    unsigned int version;
    std::string name;
    stream_.read_start_object_tag(name, version); // may throw an exc. here
    persister_base<IOS>* p = typemgr<IOS>::getPersister_xml(name);
    if(!p) {
      std::cerr << "cannot get persister for type \"" << name << "\"" << std::endl;
      abort();
    }
    p->read_members(*this);
    stream_.read_end_object_tag();
    
    //const T* add = dynamic_cast<const T*>(p->get_object_addr());
    const T* add = (const T*)(p->get_object_addr());
    if(!add) {
      std::cerr << "objio::next_object(T* obj) Object type in file does not match requirement"
		<< std::endl;
      abort();
    }
    
    
    return const_cast< T* >(add);
 }
  
  template<class T>
  void         read(CountedRef<T>& r) const {
    void* addr;
    cout << "about to read the reference tag:" << endl;
    stream_.read_start_reference_tag(addr);
    cout << "read CountedRef: addr=" << addr << endl;
    T* pt = dynamic_cast<T*>(check_address_(addr));
    if(!pt) {
      pt = read_next_object<T>();
      register_address_(addr,pt);
    }
    r = pt;
    stream_.read_end_reference_tag();
  }
  
  
  void         read(RefCount& rc) const {
    std::string name; unsigned int version;
    stream_.read_start_object_tag(name,version);
    stream_.read_end_object_tag();
  }
  
  //  template<class T>
  //  void         read(handle<T>& h) const {
  //    stream_.read(h,*this);
  //  }


  unsigned int depth() const { return 0; }

  //  virtual void read_streamer(object_structure& os) const {
  //    stream_.read_streamer(os);
  //  }

private:
  IOS stream_;
  
  mutable std::map<void const*,RefCount*> addr_;
  
  RefCount* check_address_(void const* file_addr) const {
    cout << "check_address_: addr=" << file_addr << endl;
    std::map<void const*,RefCount*>::iterator it;
    it = addr_.find(file_addr);
    if(it==addr_.end()) {
      cout << " not found..." << endl;
      return (RefCount*)0;
    }
    else {
      cout << " foud: " << it->second << endl;
      return it->second;
    }
  }
  
  void register_address_(void const* file_addr, RefCount const* addr) const {
    std::map<void const*,RefCount*>::iterator it;
    it = addr_.find(file_addr);
    if(it==addr_.end())
      addr_[(void const*)file_addr]=const_cast<RefCount*>(addr);
    else
	// should never happen. We should throw an exc. here (!)
	std::cout << "obj_input: addr already registered!" << std::endl;
  }
};


#define define_input_operator(type)                        \
template<class IOS>                                        \
obj_input<IOS> const& operator>>(obj_input<IOS> const& oo, type& v) \
{                                                          \
  oo.read(v); return oo;                                   \
}                                                          \
 
define_input_operator(int1)
define_input_operator(uint1)
define_input_operator(int2)
define_input_operator(uint2)
define_input_operator(int4)
define_input_operator(uint4)
define_input_operator(int8)
define_input_operator(uint8)
define_input_operator(float4)
define_input_operator(float8)
define_input_operator(std::string)


#define define_input_function(type)                        \
template<class IOS>                                        \
void read(obj_input<IOS> const& oo, type& v)               \
{                                                          \
  oo.read(v);                                              \
}                                                          \

define_input_function(int1)
define_input_function(uint1)
define_input_function(int2)
define_input_function(uint2)
define_input_function(int4)
define_input_function(uint4)
define_input_function(int8)
define_input_function(uint8)
define_input_function(float4)
define_input_function(float8)
define_input_function(std::string)


template<class IOS, class T>
obj_input<IOS> const& operator>>(obj_input<IOS> const& io, std::list<T>& l)
{
  io.read(l); return io;
}


template<class IOS, class T>
obj_input<IOS> const& operator>>(obj_input<IOS> const& io, std::vector<T>& v)
{
  io.read(v); return io;
}

template<class IOS, class T, class U>
obj_input<IOS> const& operator>>(obj_input<IOS> const& io, std::map<T,U>& m)
{
  io.read(m); return io;
}


template<class IOS, class T>
obj_input<IOS> const& operator>>(obj_input<IOS> const& io, CountedRef<T>& cr)
{
  io.read(cr); return io;
}

template<class IOS>
obj_input<IOS> const& operator>>(obj_input<IOS> const& io, RefCount& cr)
{
  io.read(cr); return io;
}



template<class IOS, class T>
void read(obj_input<IOS> const& oi, T* p, unsigned long sz)
{
  oi.read(p,sz);
}


template<class IOS, class T>
void read(obj_input<IOS> const& oi, persister<T,IOS>& p)
{
  oi.read(p);
}


template<class IOS>
void read(obj_input<IOS> const& oi, persister_base<IOS>* p)
{
  oi.read(p);
}


template<class IOS, class T>
void read(obj_input<IOS> const& io, std::list<T>& l)
{
  io.read(l);
}


template<class IOS, class T>
void read(obj_input<IOS> const& io, std::vector<T>& v)
{
  io.read(v);
}

template<class IOS, class T, class U>
void read(obj_input<IOS> const& io, std::map<T,U>& m)
{
  io.read(m);
}

template<class IOS, class T>
void read(obj_input<IOS> const& io, CountedRef<T>& cr)
{
  io.read(cr);
}

template<class IOS>
void read(obj_input<IOS> const& io, RefCount& cr)
{
  io.read(cr);
}


// FIXME: this may be dangerous
template<class IOS, class T>
void write(obj_output<IOS>& oo, T const& t, const char* name=0)
{
  persister<T,IOS> pp(t);
  oo.write(pp,name);
}

template<class IOS, class T>
void read(obj_input<IOS> const& io, T& t)
{
  persister<T,IOS> pp(t);
  io.read(pp);
}

#endif

