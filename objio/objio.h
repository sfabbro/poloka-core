// -*- C++ -*-
// 
// file objio.h
// 
// 
#ifndef OBJIO_H
#define OBJIO_H

#include <string>
#include <vector>
#include <list>
#include <map>

#include "toadtypes.h"



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

  //  template<class T>
  //  void   write(persister<T>& pt) {
  //    //    persister<T> pt(&t);
  //    stream_.start_object(pt.name(), pt.version());
  //    pt.write_members(*this);
  //    stream_.end_object();
  //  }

  template<class T>
  void   write(persister<T> const& p, const char* name=0) {
    stream_.start_object(p.name(), p.version());
    p.write_members(*this);
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

  // pointers are a difficult matter.
  //  template<class T>
  //  void           write(T const* t, unsigned long sz) {
  //    stream_.write(t,sz,*this);
  //  }

  template<class T>
  void           write(std::list<T> const& l, const char* name=0) {
    stream_.start_collection(l.size(), name);
    typename std::list<T>::const_iterator it;
    for(it=l.begin();it!=l.end();it++)
      *this << *it; //      write(*it);
    stream_.end_collection();
  }

  template<class T>
  void           write(std::vector<T> const& v, const char* name=0) {
    stream_.start_collection(v.size(), name);
    typename std::vector<T>::const_iterator it;
    for(it=v.begin();it!=v.end();it++)
      *this << *it; //      write(*it);
    stream_.end_collection();
  }

  template<class T, class U>
  void           write(std::map<T,U> const& m, const char* name=0) {
    stream_.start_collection(m.size(), name);
    typename std::map<T,U>::const_iterator it;
    for(it=m.begin();it!=m.end();it++)
      *this << *it; //      write(*it);
    stream_.end_collection();
  }
  
  
  template<class T, class U>
  void           write(std::pair<T,U> const& p, const char* name=0) {
    stream_.start_collection(2, name);
    *this << p.first;     //    write(p.first);
    *this << p.second;    //    write(p.second);
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
  
  template<class T> friend class persister;
};



template<class IOS, class T>
inline obj_output<IOS>& operator<<(obj_output<IOS>& oo, T t)
{
  oo.write(t);
  return oo;
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
      read(val);
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
      read(val);
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
    read(v1); p.first=v1;
    read(v2); p.second=v2;
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


template<class IOS, class T>
inline obj_input<IOS> const& operator>>(obj_input<IOS> const& oi, T& t)
{
  oi.read(t);
  return oi;
}





#endif

