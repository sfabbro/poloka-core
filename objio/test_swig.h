// -*- C++ -*-
// 
// file test_swig.h
// 
// 
#ifndef TEST_SWIG_H
#define TEST_SWIG_H


#include<string>

#include "countedref.h"
#include "persistence.h"



struct A {};
struct C {};

struct AA : public A, public C {
  CLASS_VERSION(AA,1);
  #define AA_is_persistent
};


class Point : public RefCount {
  CLASS_VERSION(Point,1);
  #define Point_is_persistent
public:
  Point() : x_(0), y_(0) {}
  virtual ~Point() {}
  
  double  x() const { return x_; }
  double  y() const { return y_; }
  double& x()       { return x_; }
  double& y()       { return y_; }

  void     print() const { std::cout << "Point: x=" << x_ << " y=" << y_ << std::endl; }  
  
protected:
  float8 x_;
  float8 y_;
};


class Star : public Point {
  CLASS_VERSION(Star,1);
  #define Star_is_persistent
public:
  Star() : Point(), id_(0), flux_(0) {}
  virtual ~Star() {}

  unsigned int  id() const { return id_; }
  unsigned int& id()       { return id_; }
  double   flux()   const { return flux_; }
  double&  flux()         { return flux_; }
  
  void     print() const { std::cout << "Star: x=" << x_ << " y=" << y_ << " flux=" << flux_ << " id=" << id_ << std::endl; }
private:
  uint4 id_;
  float8 flux_;
};


// define_template_args B<T>
// make_persister_for   B<int>
template<class T>
class B : public A {
  CLASS_VERSION(B,1);
  #define B_is_persistent
public:
  B() {}
  virtual ~B() {}
  
private:
  T t_;
};


// define_template_args BB<T,U>
// make_persister_for BB<string,string>
// make_persister_for BB<string,double>
// make_persister_for BB<double,float>
// make_persister_for BB<double,double>
template<class T, class U>
class BB : public A {
  CLASS_VERSION(BB,2);
  #define BB_is_persistent
public:
  BB() {}
  virtual ~BB() {}
  
  std::list<T>  lt_;
  std::map<  T, U> mtu_;
  
private:  
  T t_;
  U u_;
};



class Tutu {
  CLASS_VERSION(Tutu,7);
  #define Tutu_is_persistent
public:
  double toto;
  double glop[220][22][55];
  float  gloups[20];
};


#endif


