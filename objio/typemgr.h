// -*- C++ -*-
// 
// $Id: typemgr.h,v 1.4 2004/03/08 17:41:02 guy Exp $
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
  
  static persister_base<IOS>*   getPersister_rtti(std::string const& obj_rttiname) {
    typename std::map<std::string,persister_base<IOS>*>::const_iterator it;
    it = persisterMap_rtti_->find(obj_rttiname);
    if(it==persisterMap_rtti_->end()) return 0;
    return it->second;
  }
  static persister_base<IOS>*   getPersister_xml(std::string const& obj_xmlname) {
    typename std::map<std::string,persister_base<IOS>*>::const_iterator it;
    it = persisterMap_xml_->find(obj_xmlname);
    if(it==persisterMap_xml_->end()) return 0;
    return it->second;
  }
  
private:
  static std::map<std::string,persister_base<IOS>*>* persisterMap_rtti_;
  static std::map<std::string,persister_base<IOS>*>* persisterMap_xml_;
  template<class T, class U> friend class persister;
  template<class T, class U> friend class type_registrar;
};



template<class T,class IOS>
class type_registrar {
public:
  type_registrar(const std::string& classname) { register_type(classname); }
  
private:
  static void register_type(const std::string& classname) {
    std::string rttiname = typeid(T).name();
    typename std::map<std::string,persister_base<IOS>*>::const_iterator it;
    it = typemgr<IOS>::persisterMap_rtti_->find(rttiname);
    if(it==typemgr<IOS>::persisterMap_rtti_->end()) 
      (*typemgr<IOS>::persisterMap_rtti_)[rttiname] = (persister_base<IOS>*)new persister<T,IOS>;
    
    it = typemgr<IOS>::persisterMap_xml_->find(classname);
    if(it==typemgr<IOS>::persisterMap_xml_->end()) 
      (*typemgr<IOS>::persisterMap_xml_)[classname] = (persister_base<IOS>*)new persister<T,IOS>;
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
