// This may look like C code, but it is really -*- C++ -*-
#ifndef FRAME__H
#define FRAME__H

#include <iostream>
#include "point.h"

#include "../objio/persister.h"
#include "../objio/toadtypes.h"

class Gtransfo;
typedef enum {WholeSizeFrame, ClippedSizeFrame} WhichFrame;



class FitsHeader;
class Image;

//! rectangle with sides parallel to axes.
/*! when Frame's are used to define subparts of images, xMin and xMax refer
  to the first and last pixels in the subimage */
class Frame {
public:
  //! coordinate of boundary.
  float8 xMin,xMax,yMin,yMax;
  
  //! Default constructor
  Frame();
  
  //! actual image frame.
  Frame(const Image &image);
  
  //! this one is dangerous: you may swap the 2 middle arguments. 
  //! Prefer next one
  Frame(const double &xMin, const double &yMin, 
	const double &xMax, const double &yMax);
  
  //! typical use: Frame(Point(xmin,ymin),Point(xmax,ymax))
  Frame(const Point &LowerLeft, const Point &UpperRight);

  //! 2 kinds of bounds in headers, the chip size and some 
  //! that may be added by hand. See WriteInHeader().
  Frame(const FitsHeader &header, WhichFrame which=ClippedSizeFrame);
  
  
  //! number of pixels in x direction  
  double Nx() const {return xMax-xMin+1;}
  
  //! number of pixels in y direction  
  double Ny() const {return yMax-yMin+1;}
  
  //! middle of the frame
  Point Center() const {return Point((xMax+xMin)*0.5,(yMax+yMin)*0.5);}
  
  //! assumes that Transfo is a shift or involves a 'simple rotation'
  Frame ApplyTransfo(const Gtransfo &T) const;

  //! intersection of Frame's.
  Frame operator*(const Frame &Right) const;  /* intersection : a = b n c */
  
  //! intersection of Frame's
  Frame& operator*=( const Frame &Right);     /* intersection : a = a n b */
  
  //! union of Frames
  Frame operator+(const Frame &Right) const;  /* union : a = b u c */

  //! union of Frames
  Frame& operator+=( const Frame &Right);     /* intersection : a = a u b */
  
  //! shrinks the frame (if MarginSize>0), enlarges it (if MarginSize<0).
  void CutMargin(const double MarginSize);
  
  //! shrinks the frame (if MarginSize>0), enlarges it (if MarginSize<0).
  void CutMargin(const double MarginX, const double MarginY);
  
  //! necessary for comparisons (!= is defined from this one implicitely)
  bool operator ==(const Frame &Right) const;
  
  //! comparison
  bool operator !=(const Frame &Right) const {return !(*this == Right);}

  //! rescale it. The center does not move.
  Frame Rescale(const double Factor) const;
  
  // the area.
  double Area() const;
  
  //! inside?
  bool InFrame(const double &x, const double &y) const;
  
  //! same as above
  bool InFrame(const Point &pt) const 
  {return InFrame(pt.x,pt.y);}
  
  //! distance to closest boundary.
  double MinDistToEdges(const Point &P) const;
  
  //! write a frame in the fits header.
  void WriteInHeader(FitsHeader &Head) const;
  
  void dump(ostream & stream = cout) const;
  
  //! allows \verbatim cout << frame; \endverbatim.
  friend ostream & operator<<(ostream &stream, const Frame &Right) 
          { Right.dump(stream); return stream;};
  
private:
  void order();
  
  static const unsigned short __version__=2;
  friend class persister<Frame>;
};


bool HasClippedFrame(const FitsHeader &Head);

#endif /* FRAME__H */
