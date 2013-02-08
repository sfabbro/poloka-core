#ifndef FATPOINT__H
#define FATPOINT__H

#include <poloka/point.h>

//! A Point with uncertainties
class FatPoint : public Point
{
 public :
  double vx,vy,vxy;

  FatPoint() : Point() {vx=vy=1; vxy =0;};

    FatPoint(const Point &P, double Vx = 1, double Vy =1 , double Vxy=0) 
      : Point(P),  vx(Vx), vy(Vy), vxy(Vxy) {};

  FatPoint(const double X, const double Y, 
	   const double Vx=1, const double Vy=1, const double Vxy=0)
    : Point(X,Y), vx(Vx), vy(Vy), vxy(Vxy) {};


};  



#endif
