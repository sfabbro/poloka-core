// -*- C++ -*-
// 
// $Id: typemgr.cc,v 1.2 2004/03/01 22:01:44 nrl Exp $
// 
// \file typemgr.cc
// 
// 
#include "typemgr.h"
#include "xmlstream.h"


std::map<std::string,persister_base<xmlstream>*>* typemgr<xmlstream>::persisterMap_ 
   = new std::map<std::string,persister_base<xmlstream>*>();

