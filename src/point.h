// This may look like C code, but it is really -*- C++ -*-
#ifndef POINT__H
#define POINT__H
#include <iostream>
#include <cmath>

#include "persistence.h"


using namespace std;

/*! \file */

//! A point in a plane.
class Point {
  CLASS_VERSION(Point,1);
  #define Point__is__persistent
public:
  //! coordinate
   double x,y;

  //! - contructor
   Point(const double &xx=0, const double &yy=0) : x(xx), y(yy) {};

  //! -
   double Distance(const Point& Other) const { return sqrt((x-Other.x)*(x-Other.x) + (y-Other.y)*(y-Other.y));};

  //! distance squared to Other
   double Dist2(const Point& Other) const { return ((x-Other.x)*(x-Other.x) + (y-Other.y)*(y-Other.y));};

  //! Operator can be e.g. a Gtransfo.
   template<class Operator> void Apply(const Operator &Op) { *this = Op(*this);}

  //! Sum
  Point operator + (const Point &Right) const { return Point(x+Right.x, y+Right.y);}


  //! Difference
  Point operator - (const Point &Right) const { return Point(x-Right.x, y-Right.y);}

  virtual void dump(ostream& s = cout) const { s <<" x " << x << " y " << y;}  

#ifndef SWIG
  //! -
     friend ostream& operator << (ostream& stream, const Point &point) 
    { stream << " x " << point.x << " y " << point.y; return stream;}
#endif
};

#ifndef SWIG
//#ifdef __CINT__ /* rootcint  does not recognise friend inline functions */
ostream& operator << (ostream& stream, const Point &point);
//#endif
#endif


#endif /* POINT__H */
