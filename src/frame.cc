// -*- C++ -*-
// $Id: frame.cc,v 1.2 2004/04/26 13:30:39 guy Exp $
// 
// 
// 
// Last Modified: $Date: 2004/04/26 13:30:39 $
// By:            $Author: guy $
// 
#include <iostream>


#include "frame.h"
#include "fitsimage.h"
#include "gtransfo.h"


/****************** Frame class methods ***********************/
Frame::Frame(const Point &LowerLeft, const Point &UpperRight)
{
  *this = Frame(LowerLeft.x,LowerLeft.y,UpperRight.x,UpperRight.y);
}


Frame::Frame(const double &xmin, const double &ymin, 
	     const double &xmax, const double &ymax)
{
  xMin = min(xmin,xmax); xMax = max(xmin,xmax); 
  yMin = min(ymin,ymax); yMax = max(ymin,ymax);
}


Frame::Frame()
{
  xMin = xMax = yMin = yMax =0;
}


Frame::Frame(const Image &image)
{
  xMin = 0;
  yMin = 0;
  xMax = image.Nx() - 1;
  yMax = image.Ny() - 1;
}


/* clearly this assumes that the Transfo is essentially a translation or a simple rotation (n*pi/2) could probably be improved */
Frame Frame::ApplyTransfo(const Gtransfo &T, const WhichTransformed W) const
{
  // 2 opposite corners
  double xtmin1, xtmax1, ytmin1, ytmax1;
  T.apply(xMin,yMin,xtmin1,ytmin1);
  T.apply(xMax,yMax,xtmax1,ytmax1);
  Frame fr1(min(xtmin1,xtmax1), min(ytmin1,ytmax1), 
	    max(xtmin1,xtmax1), max(ytmin1,ytmax1));
  // 2 other corners
  double xtmin2, xtmax2, ytmin2, ytmax2;
  T.apply(xMin, yMax, xtmin2, ytmax2);
  T.apply(xMax, yMin, xtmax2, ytmin2);
  Frame fr2(min(xtmin2,xtmax2), min(ytmin2,ytmax2), 
	    max(xtmin2,xtmax2), max(ytmin2,ytmax2));

  if (W == SmallFrame) return fr1*fr2;
  return fr1+fr2;
}


/* positive if inside, negative if outside */
double Frame::MinDistToEdges(const Point &P) const
{
  return min(min(P.x - xMin, xMax - P.x) /* minx */, 
	     min(P.y - yMin, yMax - P.y) /* miny */);
}


bool HasClippedFrame(const FitsHeader &Head)
{
  return Head.HasKey("XNEWBEG");
}


Frame::Frame(const FitsHeader &header, WhichFrame which)
{
  if(which == ClippedSizeFrame) {
    xMin = header.KeyVal("XNEWBEG");
    yMin = header.KeyVal("YNEWBEG");
    xMax = header.KeyVal("XNEWEND");
    yMax = header.KeyVal("YNEWEND");
  }
  if(which == WholeSizeFrame ||!header.HasKey("XNEWBEG")) {
    xMin = yMin = 0;
    xMax = header.KeyVal("NAXIS1");
    xMax -=1;
    yMax = header.KeyVal("NAXIS2");
    yMax -=1;
  }
}


Frame Frame::operator*(const Frame &Right) const
{
  Frame result = *this;
  result *= Right;
  return result;
}


Frame& Frame::operator*=( const Frame &Right)
{
  Frame rightCopy = Right;
  // make sure that coordinates are properly ordered 
  this->order();
  rightCopy.order();
  xMin = max(xMin,rightCopy.xMin);
  xMax = min(xMax,Right.xMax);
  yMin = max(yMin,Right.yMin);
  yMax = min(yMax,Right.yMax);
  // check for an actual overlap. Why was this check added?
  if (xMin>xMax || yMin > yMax) *this = Frame();
  return *this;
}


Frame Frame::operator+(const Frame &Right) const
{
  Frame result = *this;
  result += Right;
  return result;
}


Frame& Frame::operator+=( const Frame &Right)
{
  Frame rightCopy = Right;
  // make sure that coordinates are properly ordered 
  this->order();
  rightCopy.order();
  xMin = min(xMin, rightCopy.xMin);
  xMax = max(xMax, Right.xMax);
  yMin = min(yMin, Right.yMin);
  yMax = max(yMax, Right.yMax);
  return *this;
}

void Frame::order()
{
  if (xMin > xMax ) swap (xMin,xMax);
  if (yMin > yMax ) swap (yMin,yMax);
}


bool  Frame::operator==( const Frame &Right) const
{
  return ((xMin == Right.xMin) && (xMax == Right.xMax) &&
	  (yMin == Right.yMin) && (yMax == Right.yMax));
}


void Frame::CutMargin(const double MarginX, const double MarginY)
{
  xMin += MarginX;
  yMin += MarginY;
  xMax -= MarginX;
  yMax -= MarginY;
}


void Frame::CutMargin(const double MarginSize)
{
  CutMargin(MarginSize, MarginSize);
}


Frame Frame::Rescale(const double Factor) const
{
  double hxsize = fabs(Factor*0.5*(xMax - xMin));
  double xcenter = 0.5*(xMax + xMin);
  double hysize = fabs(Factor*0.5*(yMax - yMin));
  double ycenter = 0.5*(yMax + yMin);
  return Frame(xcenter - hxsize , ycenter - hysize, 
	       xcenter + hxsize, ycenter + hysize);
}


double Frame::Area() const
{
  return fabs((xMax - xMin)*(yMax - yMin));
}


bool Frame::InFrame(const double &x,const double &y) const
{
  return ((x <= xMax) && (y<=yMax) && 
	  (x>=xMin) && (y>=yMin));
}


void Frame::WriteInHeader(FitsHeader &Head) const /* could be in the FitsHeader class as well */
{
  Head.AddOrModKey("XNEWBEG",xMin,
		 "new Xbeg of image after geometric transfo.");
  Head.AddOrModKey("YNEWBEG",yMin,
		 "new Ybeg of image after geometric transfo.");
  Head.AddOrModKey("XNEWEND",xMax,
		 "new Xend of image after geometric transfo.");
  Head.AddOrModKey("YNEWEND",yMax,
		 "new Yend of image after geometric transfo.");
}


void  Frame::dump(ostream & stream) const
{
  stream << "xmin ymin "  << xMin << ' ' << yMin 
	 << " xmax ymax " << xMax << ' ' << yMax << endl;
}

