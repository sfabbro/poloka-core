// -*- C++ -*-
// 
// $Id: objio.h,v 1.4 2004/02/27 16:34:36 nrl Exp $
// 
// \file objio.h
// 
// 
#ifndef OBJIO_H
#define OBJIO_H

#include <string>
#include <vector>
#include <list>
#include <map>

#include "toadtypes.h"


class persister_base;
template<class T> class persister;


template<class IOS>
class obj_output { //: public obj_output_base {
 public:
  obj_output(const std::string& name) : stream_(name) {}
  ~obj_output() { }

  void         open(std::string const& filename, int compression) const {
    stream_.open(filename,compression);
  }

  void         close() { stream_.close(); }

  template<class T>
  void   write(persister<T> const& p, const char* name=0) {
    stream_.start_object(p.name(), p.version());
    p.write_members(*this);
    stream_.end_object();
  }
  
  void   write(persister_base const* p, const char* name=0) {
    stream_.start_object(p->name(), p->version());
    p->write_members(*this);
    stream_.end_object();
  }
  
  virtual void   write(int1 v,   const char* name=0) { stream_.write(v, name); }
  virtual void   write(uint1 v,  const char* name=0) { stream_.write(v, name); }
  virtual void   write(int2 v,   const char* name=0) { stream_.write(v, name); }
  virtual void   write(uint2 v,  const char* name=0) { stream_.write(v, name); }
  virtual void   write(int4 v,   const char* name=0) { stream_.write(v, name); }
  virtual void   write(uint4 v,  const char* name=0) { stream_.write(v, name); }
  virtual void   write(int8 v,   const char* name=0) { stream_.write(v, name); }
  virtual void   write(uint8 v,  const char* name=0) { stream_.write(v, name); }
  virtual void   write(float4 v, const char* name=0) { stream_.write(v, name); }
  virtual void   write(float8 v, const char* name=0) { stream_.write(v, name); }
  virtual void   write(const std::string& v, const char* name=0) { stream_.write(v, name); }
  
  
  // tabs     //    stream_.write(t,sz,*this);
  /*
    template<class T>
    void           write(T const* t, unsigned long sz, const char* name=0) {
    stream_.start_collection(name,sz,t);
    int i;
    for(i=0;i<sz;i++) write(t+i)
    stream_.end_collection();
    }
  */
  
  /*
  // raw pointers 
  template<class T>
  void           write(T const* t) {
    bool w = true; //check_pointer_((void*)t); // did we write this object already ?
    persister_base* b = objmgr::getPersister(typeinfo(t).name());
    stream_.start_pointer(name);
    if(!w) {
      *b = t;
      //      *this << *b;
    }
    else 
      stream_.write(t);
    stream_.end_pointer();
  }
  */
  
  template<class T>
  void           write(std::list<T> const& l, const char* name=0) {
    stream_.start_collection(l.size(), name);
    typename std::list<T>::const_iterator it;
    for(it=l.begin();it!=l.end();it++)
      *this << *it;
    stream_.end_collection();
  }

  template<class T>
  void           write(std::vector<T> const& v, const char* name=0) {
    stream_.start_collection(v.size(), name);
    typename std::vector<T>::const_iterator it;
    for(it=v.begin();it!=v.end();it++) {
      *this << *it;
    }
    stream_.end_collection();
  }

  template<class T, class U>
  void           write(std::map<T,U> const& m, const char* name=0) {
    stream_.start_collection(m.size(), name);
    typename std::map<T,U>::const_iterator it;
    for(it=m.begin();it!=m.end();it++)
      // this time, we must use write().
      // Because we know what we are writing out...
      write(*it);
    stream_.end_collection();
  }
  
  
  template<class T, class U>
  void           write(std::pair<T,U> const& p, const char* name=0) {
    stream_.start_collection(2, name);
    *this << p.first;
    *this << p.second;
    stream_.end_collection();
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
  //  dict_output<IOS>& d_output_;
  
  bool check_pointer(void*);
  
  template<class T> friend class persister;
};


#define define_output_operator(type)                       \
template<class IOS>                                        \
obj_output<IOS>& operator<<(obj_output<IOS>& oo, type v)   \
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


template<class IOS>
class obj_input { //: public obj_input_base {
public:
  obj_input(const std::string& name) : stream_(name) {}
  ~obj_input() {}
  
  void         open(std::string const& filename) {
    stream_.open(filename);
  }
  
  void         close() { stream_.close(); }
  
  //  template<class T>
  //  virtual persister_base*  read() const;

  //  virtual persister_base* read() const { return 0; }

  template<class T>
  void         read(persister<T>& p) const {
    unsigned int version;
    std::string name;
    stream_.start_object(name, version); // may throw an exc. here
    p.read_members(*this);
    stream_.end_object();
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

  //  template<class T>
  //  void         read(T*& t, unsigned long& sz) const {
  //    stream_.read(t,sz,*this);
  //  }

  template<class T>
  void         read(std::list<T>& l) const {
    unsigned int i, sz;
    T val;
    stream_.start_collection(sz);
    for(i=0;i<sz;i++) {
      //      read(val);
      *this >> val;
      l.push_back(val);
    }
    stream_.end_collection();
  }

  template<class T>
  void         read(std::vector<T>& v) const {
    unsigned int i, sz;
    T val;
    stream_.start_collection(sz);
    for(i=0;i<sz;i++) {
      //      read(val);
      *this >> val;
      v.push_back(val);
    }
    stream_.end_collection();
  }
  
  template<class T, class U>
  void         read(std::map<T,U>& m) const {
    unsigned int i, sz;
    std::pair<T,U> val;
    stream_.start_collection(sz);
    for(i=0;i<sz;i++) {
      read(val);
      m[val.first]=val.second;
    }
    stream_.end_collection();
  }
  
  template<class T, class U>
  void         read(std::pair<T,U>& p) const {
    unsigned int sz;
    T v1; U v2;
    stream_.start_collection(sz);
    *this >> v1; p.first=v1;
    *this >> v2;
    p.second=v2;
    stream_.end_collection();
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




#endif

