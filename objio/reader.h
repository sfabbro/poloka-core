// -*- C++ -*-
// 
// file reader.h
// 
// complex and simple type readers.
//
#ifndef READER_H
#define READER_H

#include <list>

#include "objio.h"
#include "dict.h"
#include "persister.h"


template<class O, class OIS>
class reader_base {
public:
  virtual void read(obj_input<OIS> const&,O&) {}
};


template<class O, class OIS>
class class_reader : public reader_base<O, OIS> {
public:
  class_reader() { }
  ~class_reader() { clear_(); }
  
  virtual void read(obj_input<OIS> const& oi, O& o) {
    typename std::list<reader_base<O,OIS>*>::iterator it;
    for(it=rdlist_.begin();it!=rdlist_.end();it++)
      it->read(oi, o);
  }
  
  void         link(dict const& d, persister<O> const& p);
  
private:
  std::list<reader_base<O,OIS>*> rdlist_;
  
  void clear_() { rdlist_.clear(); }
};


template<class O, class OIS, class T, class D>
class simple_type_reader : public reader_base<O, OIS> {
public:
  simple_type_reader(D O::* p) : p_(p) { }
  ~simple_type_reader() { }
  
  virtual void read(obj_input<OIS> const& oi, O& o) {
    T val; 
    oi.read(val);
    o.*p_ = (D)val;
  }

private:
  D O::* p_;
};


template<class O, class OIS, class T>
class complex_type_reader : public reader_base<O, OIS> {
public:
  complex_type_reader(T O::* p) : p_(p) { }
  ~complex_type_reader() { }
  
  virtual void read(obj_input<OIS> const& oi, O& o) { oi.read(o.*p_); }

private:
  T O::* p_;
};



//////////////////////////// INLINE METHODS ///////////////////////////
template<class O, class OIS>
void  class_reader<O,OIS>::link(dict const& d, persister<O> const& p)
{
  clear_();
  
  unsigned int i, class_i;
  std::string dict_memberName;
  std::string dict_memberType;
  std::string class_memberType;
  
  for(i=0;i<d.size();i++) {
    dict_memberName = d.member(i).name;
    dict_memberType = d.member(i).typeName;
    
    if(!p.hasName(memberName)) {
      //      rdlist_.push_back(ReaderMgr::getSkip(dict_memberType));
      continue;
    }
    class_i = p.getNameIndex(memberName);
    class_memberType = p.memberType(class_i);
    //    rdlist_.push_back(ReaderMgr::getReader(class_memberType, class_memberName));
  }
}

#endif

