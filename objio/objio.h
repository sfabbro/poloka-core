// -*- C++ -*-
// $Id: objio.h,v 1.6 2004/03/02 12:40:16 nrl Exp $
// 
// \file objio.h
// 
// Last modified: $Date: 2004/03/02 12:40:16 $
// by:            $Author: nrl $
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
    check_address_(p.obj_);
    stream_.write_start_object_tag(p.name(), p.version(), p.get_object_addr());
    p.write_members(*this);
    stream_.write_end_object_tag();
  }
  
  void   write(persister_base<IOS> const* p, const char* name=0) {
    check_address_(p->get_object_addr());
    stream_.write_start_object_tag(p->name(), p->version(), p->get_object_addr());
    p->write_members(*this);
    stream_.write_end_object_tag();
  }
  
  virtual void   write(int1 const& v,   const char* name=0) { stream_.write(v, name, &v); }
  virtual void   write(uint1 const& v,  const char* name=0) { stream_.write(v, name, &v); }
  virtual void   write(int2 const& v,   const char* name=0) { stream_.write(v, name, &v); }
  virtual void   write(uint2 const& v,  const char* name=0) { stream_.write(v, name, &v); }
  virtual void   write(int4 const& v,   const char* name=0) { stream_.write(v, name, &v); }
  virtual void   write(uint4 const& v,  const char* name=0) { stream_.write(v, name, &v); }
  virtual void   write(int8 const& v,   const char* name=0) { stream_.write(v, name, &v); }
  virtual void   write(uint8 const& v,  const char* name=0) { stream_.write(v, name, &v); }
  virtual void   write(float4 const& v, const char* name=0) { stream_.write(v, name, &v); }
  virtual void   write(float8 const& v, const char* name=0) { stream_.write(v, name, &v); }
  virtual void   write(const std::string& v, const char* name=0) { stream_.write(v, name, &v); }
  
  
  // tabs     //    stream_.write(t,sz,*this);
  template<class T>
  void           write(T const* t, unsigned long sz, const char* name=0) {
    stream_.write_start_collection_tag(sz,name,t);
    int i;
    for(i=0;i<sz;i++)
      write(*(t+i));
    stream_.write_end_collection_tag();
  }
  
  // raw pointers 
  template<class T>
  void           write(T const* t, const char* name=0) {
    bool w = check_address_((void const*)t); // did we write this object already ?
    persister<T,IOS>* b = (persister<T,IOS>*)typemgr<IOS>::getPersister(typeid(*t).name());
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
    stream_.write_start_collection_tag(l.size(), name,&l);
    typename std::list<T>::const_iterator it;
    for(it=l.begin();it!=l.end();it++)
      *this << *it;
    stream_.write_end_collection_tag();
  }

  template<class T>
  void           write(std::vector<T> const& v, const char* name=0) {
    stream_.write_start_collection_tag(v.size(), name,&v);
    typename std::vector<T>::const_iterator it;
    for(it=v.begin();it!=v.end();it++) {
      *this << *it;
    }
    stream_.write_end_collection_tag();
  }

  template<class T, class U>
  void           write(std::map<T,U> const& m, const char* name=0) {
    stream_.write_start_collection_tag(m.size(), name,&m);
    typename std::map<T,U>::const_iterator it;
    for(it=m.begin();it!=m.end();it++)
      // this time, we must use write().
      // Because we know what we are writing out...
      write(*it);
    stream_.write_end_collection_tag();
  }
  
  
  template<class T, class U>
  void           write(std::pair<T,U> const& p, const char* name=0) {
    stream_.write_start_collection_tag(2, name,&p);
    *this << p.first;
    *this << p.second;
    stream_.write_end_collection_tag();
  }
  
  //  template<class T>
  //  void         write(CountedRef<T> const&) {}

  //  template<class T>
  //    void           write(handle<T> const& h) {
  //    stream_.write(h,*this);
  //  }

  //  virtual void   write_streamer(streamer_base const* s) {
  //    stream_.write_streamer(s);
  //  }
  
  // 
  // dictionaries
  // 
  

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
    for(i=0;i<sz;i++) {
      //      read(val);
      *this >> val;
      l.push_back(val);
    }
    stream_.read_end_collection_tag();
  }

  template<class T>
  void         read(std::vector<T>& v) const {
    unsigned int i, sz;
    T val;
    stream_.read_start_collection_tag(sz);
    for(i=0;i<sz;i++) {
      //      read(val);
      *this >> val;
      v.push_back(val);
    }
    stream_.read_end_collection_tag();
  }
  
  template<class T, class U>
  void         read(std::map<T,U>& m) const {
    unsigned int i, sz;
    std::pair<T,U> val;
    stream_.read_start_collection_tag(sz);
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
    *this >> v1; p.first=v1;
    *this >> v2;
    p.second=v2;
    stream_.read_end_collection_tag();
  }
  
  //  template<class T>
  //  void         read(CountedRef<T>&) const {}

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



#endif

