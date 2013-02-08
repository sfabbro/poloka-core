// This may look like C code, but it is really -*- C++ -*-
#ifndef POINT__H
#define POINT__H
#include <iostream>
#include <cmath>


using namespace std;

/*! \file */

//! A point in a plane.
class Point {
public:
  virtual ~Point() {}

  //! coordinate
   double x,y;

  //! - contructor
   Point() : x(0), y(0) {};

  //! - contructor
   Point(const double &xx, const double &yy) : x(xx), y(yy) {};

  //! -
   double Distance(const Point& Other) const { return sqrt((x-Other.x)*(x-Other.x) + (y-Other.y)*(y-Other.y));};

  //! distance squared to Other
   double Dist2(const Point& Other) const { return ((x-Other.x)*(x-Other.x) + (y-Other.y)*(y-Other.y));};

  //! Sum
  Point operator + (const Point &Right) const { return Point(x+Right.x, y+Right.y);}


  //! Difference
  Point operator - (const Point &Right) const { return Point(x-Right.x, y-Right.y);}

  virtual void dump(ostream& s = cout) const { s <<" x " << x << " y " << y;}  

  //! -
     friend ostream& operator << (ostream& stream, const Point &point) 
    { stream << " x " << point.x << " y " << point.y; return stream;}
};

ostream& operator << (ostream& stream, const Point &point);

#endif /* POINT__H */
