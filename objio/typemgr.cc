// -*- C++ -*-
// 
// $Id: typemgr.cc,v 1.3 2004/03/08 17:41:02 guy Exp $
// 
// \file typemgr.cc
// 
// 
#include "typemgr.h"
#include "xmlstream.h"


std::map<std::string,persister_base<xmlstream>*>* typemgr<xmlstream>::persisterMap_rtti_ 
   = new std::map<std::string,persister_base<xmlstream>*>();

std::map<std::string,persister_base<xmlstream>*>* typemgr<xmlstream>::persisterMap_xml_ 
   = new std::map<std::string,persister_base<xmlstream>*>();


