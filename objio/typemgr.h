// -*- C++ -*-
// 
// $Id: typemgr.h,v 1.1 2004/02/27 16:30:35 nrl Exp $
// 
// \file typemgr.h
// 
// 
#ifndef TYPEMGR_H
#define TYPEMGR_H

#include <typeinfo>
#include <string>

#include "persister.h"

// typemgr maintains a static map<> of persister rttiName 
// and persister creators. Each persister is supposed
// to register itself with typemgr 
class typemgr {
public:
  
  static int    size() { return persisterMap_->size(); }
  
  static persister_base*     getPersister(std::string const& obj_pointer_rttiname) {
    std::map<std::string,persister_base*>::const_iterator it;
    it = persisterMap_->find(obj_pointer_rttiname);
    if(it==persisterMap_->end()) return 0;
    return it->second;
  }
  
private:
  static std::map<std::string,persister_base*>* persisterMap_;
  
  template<class T> friend class persister;
  template<class T> friend class type_registrar;
};



template<class T>
class type_registrar {
public:
  type_registrar() { register_type(); }
  
private:
  static void register_type() {
    std::string rttiname = typeid(T*).name();
    std::map<std::string,persister_base*>::const_iterator it;
    it = typemgr::persisterMap_->find(rttiname);
    if(it==typemgr::persisterMap_->end()) 
      (*typemgr::persisterMap_)[rttiname] = (persister_base*)new persister<T>;
  }
};


#endif

