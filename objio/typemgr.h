// -*- C++ -*-
// 
// $Id: typemgr.h,v 1.2 2004/03/01 22:01:44 nrl Exp $
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
template<class IOS>
class typemgr {
public:
  
  static int    size() { return persisterMap_->size(); }
  
  static persister_base<IOS>*   getPersister(std::string const& obj_rttiname) {
    std::map<std::string,persister_base<IOS>*>::const_iterator it;
    it = persisterMap_->find(obj_rttiname);
    if(it==persisterMap_->end()) return 0;
    return it->second;
  }
  
private:
  static std::map<std::string,persister_base<IOS>*>* persisterMap_;
  
  template<class T, class U> friend class persister;
  template<class T, class U> friend class type_registrar;
};



template<class T,class IOS>
class type_registrar {
public:
  type_registrar() { register_type(); }
  
private:
  static void register_type() {
    std::string rttiname = typeid(T).name();
    std::map<std::string,persister_base<IOS>*>::const_iterator it;
    it = typemgr<IOS>::persisterMap_->find(rttiname);
    if(it==typemgr<IOS>::persisterMap_->end()) 
      (*typemgr<IOS>::persisterMap_)[rttiname] = (persister_base<IOS>*)new persister<T,IOS>;
  }
};


#endif



#ifdef GARBAGE
  static persister_base<IOS>*   getPersister(void* t) {
    std::map<std::string,persister_base<IOS>*>::const_iterator it;
    std::string obj_rttiname = typeid(*t).name();
    std::cout << obj_rttiname << std::endl;
    it = persisterMap_->find(obj_rttiname);
    if(it==persisterMap_->end()) return 0;
    return it->second->clone(t);
  }
#endif
