// -*- C++ -*-
// $Id: persister.h,v 1.6 2004/04/15 11:42:17 guy Exp $
// 
// \file persister.h
// A persister is a handle on an object, which knows
// how to read this object from an input stream, or
// how to write it to an output stream. The stream 
// is generic and can be implemented many ways.
// 
// Last modified: $Date: 2004/04/15 11:42:17 $
// by:            $Author: guy $
// 
#ifndef PERSISTER_H
#define PERSISTER_H

#include "objio.h"


class obj_intput_base;
class obj_output_base;


template<class IOS>
class persister_base {
public:
  
  virtual unsigned int     version() const=0;
  virtual std::string      name() const=0;
  virtual persister_base*  clone() const=0;
  virtual void const*      get_object_addr() const=0;
  
protected:
  virtual void             write_members(obj_output<IOS>& oo) const {}
  virtual void             read_members(obj_input<IOS> const& oi) {}
  
  template<class U> friend class obj_input;
  template<class U> friend class obj_output;
};



template<class T>
class handle {
public:
  handle() : obj_(new T), own_obj_(true) { }
  handle(T* t) : obj_(t), own_obj_(false) { }
  handle(handle<T> const& h) : obj_(0) { copy(h); }
  ~handle() { if(obj_&&own_obj_) delete obj_; };
  
  void               set_object(T const* t) {  // we really need const_handles ...
    clear_(); 
    obj_=const_cast<T*>(t); own_obj_=false; 
  }
  
  inline handle<T>&  operator=(T* t) {
    clear_(); 
    obj_=t; own_obj_=false; 
    return *this; 
  }
  
  inline T*          operator->() { return obj_; }
  inline T&          operator*() { return *obj_; }
  inline T const&    operator*() const { return *obj_; }
  
  inline void        copy(handle<T> const& h) { 
    obj_=h.obj_; 
    own_obj_=false; 
  }
  
  inline handle<T> const& operator=(handle<T> const& t) { 
    copy(t); 
    return *this; 
  }

protected:
  T*   obj_;
  bool own_obj_;
  
  inline void clear_() { 
    if(obj_&&own_obj_) delete obj_; 
    obj_=0; 
    own_obj_=false; 
  }
};



// Default persister (empty)
template<class T, class IOS>
class persister : public handle<T> {
public:
  persister() : handle<T>() {}
  persister(T& obj) : handle<T>(&obj) {}
  persister(T const & obj) : handle<T>(const_cast<T*>(&obj)) {}// I know, that's UGLY.
  persister(persister<T,IOS> const& p) : handle<T>(p.obj_) {}
  virtual ~persister() {}
  
  unsigned int         version() const { return 0; }
  std::string          name() const { return (std::string)"non persistent member"; }
  unsigned int         size() const { return 0; }
  std::string          name(unsigned int i) const { return (std::string)""; }
  std::string          type(unsigned int i) const { return (std::string)""; }
  virtual void const*  get_object_addr() const { return (void*)obj_; }
  
protected:
  virtual void         write_members(obj_output<IOS>& oo) const {}
  virtual void         read_members(obj_input<IOS> const& oi) {}

  template<class U> friend class obj_input;
  template<class U> friend class obj_output;
};


#endif
