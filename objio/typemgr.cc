// -*- C++ -*-
// 
// $Id: typemgr.cc,v 1.1 2004/02/27 16:30:35 nrl Exp $
// 
// \file typemgr.cc
// 
// 
#include "typemgr.h"


std::map<std::string,persister_base*>* typemgr::persisterMap_ = new std::map<std::string,persister_base*>();

