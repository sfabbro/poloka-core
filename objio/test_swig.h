// -*- C++ -*-
// 
// file test_swig.h
// 
// 
#ifndef TEST_SWIG_H
#define TEST_SWIG_H


#include<string>
#include "persister.h"
#include "objio.h"


#define CLASS_VERSION(className,id) \
 static const unsigned short __version__=id; \
 friend class persister<className>;


struct A {};
struct C {};
struct AA : public A, public C {
private:
CLASS_VERSION(AA,1);
//static const unsigned short __version__=2;
};


class Point {
public:
  Point() : x_(0), y_(0) {}
  ~Point() {}
  
  double  x() const { return x_; }
  double  y() const { return y_; }
  double& x()       { return x_; }
  double& y()       { return y_; }

  void     print() const { std::cout << "Point: x=" << x_ << " y=" << y_ << std::endl; }  
  
protected:
  float8 x_;
  float8 y_;
  CLASS_VERSION(Point,1);
};


class Star : public Point {
public:
  Star() : Point(), id_(0), flux_(0) {}
  ~Star() {}

  unsigned int  id() const { return id_; }
  unsigned int& id()       { return id_; }
  double   flux()   const { return flux_; }
  double&  flux()         { return flux_; }
  
  void     print() const { std::cout << "Star: x=" << x_ << " y=" << y_ << " flux=" << flux_ << " id=" << id_ << std::endl; }
private:
  uint4 id_;
  float8 flux_;
  CLASS_VERSION(Star,1);
};


//template<class T>
//class B : public A {
//public:
//  B() {}
//  ~B() {}
//  
//private:
//  T t_;
//};


#endif


