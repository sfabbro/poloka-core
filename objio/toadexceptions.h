// -*- C++ -*-
// 
// 
// file toadexceptions.h 
// 
/*!
  \file toadexceptions.h
  \brief Base classes and macros for toad exceptions
    
 */
#ifndef TOADEXCEPTIONS_H
#define TOADEXCEPTIONS_H

using namespace std;
#include <string>
#include <sstream>

//namespace toads {

#define BuildExcMsg(msg) Throwable::buildMessage(msg, __FILE__, __LINE__)



class Throwable {
public:
  //! constructor
  Throwable(string const& msg) { msg_=msg; } 
  
  //! destructor -- does nothing
  virtual ~Throwable() { } 
  
  //! append stuff to the initial message
  void   append(string const& msg) { msg_+=msg; }
  
  //! return the message
  string message() const { return msg_; }
  
  //! uniform message format
  static string buildMessage(string msg, char const* file, unsigned int line) {
    stringstream ret;
    ret << "[" << file << "]" << "{" << line << "} " << msg;
    return ret.str();
  }

private:
  string msg_;
};


////////////////////////////////////////////////////////////////////////////////

class InconsistentException : public Throwable {
public:
  //! constructor 
  InconsistentException(string const& msg) : Throwable(msg) { } 
  
  //! destructor
  ~InconsistentException() { }
};

//} //namespace 

#endif


