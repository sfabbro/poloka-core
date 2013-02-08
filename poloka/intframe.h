// This may look like C code, but it is really -*- C++ -*-
#ifndef INTFRAME__H
#define INTFRAME__H

#include <iostream>



typedef enum {WholeSizeIntFrame, ClippedSizeIntFrame} WhichIntFrame;


//typedef enum {LargeIntFrame, SmallIntFrame} WhichTransformed;
class FitsHeader;
class Point;
class IntPoint;

//! rectangle with sides parallel to axes.
/*! used to hold limits of rectangles. When used to label pixels in images,
  pixels have xmin<=x<xmax and the same for y */
class IntFrame {
  
public:

  int  xmin,xmax,ymin,ymax, nx ;
  
  //! Default constructor
  IntFrame();
  
  
  //! this one is dangerous: you may swap the 2 middle arguments. 
  //! Prefer next one
  IntFrame(const int &xmin, const int &ymin, 
	   const int &xmax, const int &ymax);
  
  //! 2 kinds of bounds in headers, the chip size and some 
  //! that may be added by hand. See WriteInHeader().
  IntFrame(const FitsHeader &header, WhichIntFrame which=ClippedSizeIntFrame);
  
  
  //! number of pixels in x direction  
  int Nx() const {return xmax-xmin;}
  
  //! number of pixels in y direction  
  int Ny() const {return ymax-ymin;}
  
  //! intersection of IntFrame's.
  IntFrame operator*(const IntFrame &Right) const;  /* intersection : a = b n c */
  
  //! intersection of IntFrame's
  IntFrame& operator*=( const IntFrame &Right);     /* intersection : a = a n b */
  
  //! union of IntFrames
  IntFrame operator+(const IntFrame &Right) const;  /* union : a = b u c */

  //! union of IntFrames
  IntFrame& operator+=( const IntFrame &Right);     /* intersection : a = a u b */
  
  IntFrame Shift(const int DeltaX, const int DeltaY) const;

  IntFrame Shift(const IntPoint &P) const;

  //! shrinks the intframe (if MarginSize>0), enlarges it (if MarginSize<0).
  void CutMargin(const int MarginSize);
  
  //! shrinks the intframe (if MarginSize>0), enlarges it (if MarginSize<0).
  void CutMargin(const int MarginX, const int MarginY);
  

  bool SameSize(const IntFrame &Right) const
  {
    return ((*this) == Right);
  }

  //! necessary for comparisons (!= is defined from this one implicitely)
  bool operator ==(const IntFrame &Right) const
  {
    return ((xmin == Right.xmin) && (xmax == Right.xmax) &&
	    (ymin == Right.ymin) && (ymax == Right.ymax));
  }

  
  //! comparison
  bool operator !=(const IntFrame &Right) const {return !(*this == Right);}

  // the area.
  int Ntot() const { return (xmax - xmin)*(ymax - ymin);}
  
  //! inside?
  bool IsInside(const int x, const int y) const
  {
   return ((x < xmax) && (y<ymax) && 
	  (x>=xmin) && (y>=ymin));
  }

  bool IsOutside(const int x, const int y) const
  {
    return ((x<xmin) || (x>= xmax) || (y<ymin) || (y>=ymax));
  }


  
  //! distance to closest boundary.
  double MinDistToEdges(const Point &P) const;
  
  //! write a intframe in the fits header.
  void WriteInHeader(FitsHeader &Head) const;
  
  void dump(std::ostream & stream = std::cout) const
  {
  stream << "[(" << xmin << ',' << ymin 
	 << "):(" << xmax << ',' << ymax << ")[ ";
  }

  //! allows \verbatim cout << intframe; \endverbatim.
  friend std::ostream & operator<<(std::ostream &stream, const IntFrame &Right) 
          { Right.dump(stream); return stream;};
  
private:
  void order();
};



#endif /* INTFRAME__H */
