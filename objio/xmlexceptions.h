// -*- C++ -*-
// 
// file xmlexceptions.h
// 
// 
#ifndef XMLEXCEPTIONS_H
#define XMLEXCEPTIONS_H

#include "toadexceptions.h"

class XMLException : public Throwable {
public:
  //! Default constructor
  XMLException(string const& msg) : Throwable(msg) { } 
  
  //! Destructor 
  virtual ~XMLException() { } 
};


#endif

