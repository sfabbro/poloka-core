// -*- C++ -*-
// $Id: intframe.cc,v 1.2 2006/12/22 13:35:40 guy Exp $
// 
// 
// 
// Last Modified: $Date: 2006/12/22 13:35:40 $
// By:            $Author: guy $
// 
#include <iostream>


#include "intframe.h"
#include "fitsimage.h"
#include "intpoint.h"


/****************** IntFrame class methods ***********************/

IntFrame::IntFrame(const int &Xmin, const int &Ymin, 
	     const int &Xmax, const int &Ymax)
{
  xmin = min(Xmin,Xmax); xmax = max(Xmin,Xmax); 
  ymin = min(Ymin,Ymax); ymax = max(Ymin,Ymax);
  nx = xmax-xmin;
}


IntFrame::IntFrame()
{
  xmin = xmax = ymin = ymax =0;
  nx = 0;
}


#ifdef STORAGE
/* positive if inside, negative if outside */
double IntFrame::MinDistToEdges(const Point &P) const
{
  return min(min(P.x - xmin, xmax - P.x) /* minx */, 
	     min(P.y - ymin, ymax - P.y) /* miny */);
}
#endif


IntFrame::IntFrame(const FitsHeader &header, WhichIntFrame which)
{
  if(which == ClippedSizeIntFrame) {
    xmin = header.KeyVal("XNEWBEG");
    ymin = header.KeyVal("YNEWBEG");
    xmax = header.KeyVal("XNEWEND");
    ymax = header.KeyVal("YNEWEND");
    nx = xmax-xmin;
  }
  if(which == WholeSizeIntFrame ||!header.HasKey("XNEWBEG")) {
    xmin = ymin = 0;
    xmax = header.KeyVal("NAXIS1");
    xmax -=1;
    ymax = header.KeyVal("NAXIS2");
    ymax -=1;
    nx = xmax-xmin;
  }
}


IntFrame IntFrame::operator*(const IntFrame &Right) const
{
  IntFrame result = *this;
  result *= Right;
  return result;
}


IntFrame& IntFrame::operator*=( const IntFrame &Right)
{
  IntFrame rightCopy = Right;
  // make sure that coordinates are properly ordered 
  this->order();
  rightCopy.order();
  xmin = max(xmin,rightCopy.xmin);
  xmax = min(xmax,Right.xmax);
  ymin = max(ymin,Right.ymin);
  ymax = min(ymax,Right.ymax);
  // check for an actual overlap. Why was this check added?
  if (xmin>xmax || ymin > ymax) *this = IntFrame();
  nx = xmax - xmin;
  return *this;
}


IntFrame IntFrame::operator+(const IntFrame &Right) const
{
  IntFrame result = *this;
  result += Right;
  return result;
}


IntFrame& IntFrame::operator+=( const IntFrame &Right)
{
  IntFrame rightCopy = Right;
  // make sure that coordinates are properly ordered 
  this->order();
  rightCopy.order();
  xmin = min(xmin, rightCopy.xmin);
  xmax = max(xmax, Right.xmax);
  ymin = min(ymin, Right.ymin);
  ymax = max(ymax, Right.ymax);
  nx = xmax - xmin;
  return *this;
}


IntFrame IntFrame::Shift(const int DeltaX, const int DeltaY) const
{
  return IntFrame(xmin+DeltaX, ymin + DeltaY, xmax+DeltaX, ymax+DeltaY);
}

IntFrame IntFrame::Shift(const IntPoint &Delta) const
{
  return IntFrame(xmin+Delta.x, ymin + Delta.y, xmax+Delta.x, ymax+Delta.y);
}


void IntFrame::order()
{
  if (xmin > xmax ) swap (xmin,xmax);
  if (ymin > ymax ) swap (ymin,ymax);
  nx = xmax - xmin;
}


void IntFrame::CutMargin(const int MarginX, const int MarginY)
{
  xmin += MarginX;
  ymin += MarginY;
  xmax -= MarginX;
  ymax -= MarginY;
  nx = xmax - xmin;
}


void IntFrame::CutMargin(const int MarginSize)
{
  CutMargin(MarginSize, MarginSize);
}



