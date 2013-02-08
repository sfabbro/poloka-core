// -*- C++ -*-
// $Id: frame.cc,v 1.2 2006/12/22 13:35:41 guy Exp $
// 
// 
// 
// Last Modified: $Date: 2006/12/22 13:35:41 $
// By:            $Author: guy $
// 
#include <iostream>

#include <poloka/frame.h>
#include <poloka/fitsimage.h>


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

void Frame::read(istream& stream) 
{
  string dummy;
  stream >> dummy >> dummy >> xMin >> yMin
	 >> dummy >> dummy >> xMax >> yMax;
}
