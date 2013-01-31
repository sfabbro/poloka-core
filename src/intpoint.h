#ifndef INTPOINT__H
#define INTPOINT__H


#include "point.h"

//! 2d point with integer coordinates
struct IntPoint
{
  int x,y;
  IntPoint() : x(0), y(0) {};
  IntPoint(const int X, const int Y) : x(X), y(Y) {};


  operator Point() const { return Point(double(x), double(y));};

  
  Point operator +(const Point &Right) const 
  { return Point(x+Right.x, y+Right.y);}

};


#endif /* INTPOINT__H */
