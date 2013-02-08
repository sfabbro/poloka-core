// This may look like C code, but it is really -*- C++ -*-
#ifndef DIMAGE__H
#define DIMAGE__H

#include <list>
#include <string>

#include <poloka/basestar.h>

typedef double DPixel;
class Image;
class Kernel;

//! yet another rectangle/frame/window/subimage class. used in Stamp and Vignet
struct Window { 
  Window() : xstart(0), ystart(0), xend(0), yend(0) {}
  Window(const int Xstart, const int Ystart, 
	 const int Xend, const int Yend) : xstart(Xstart), ystart(Ystart), xend(Xend), yend(Yend) {}
  int xstart, ystart, xend, yend; 
  int Nx() const { return xend-xstart; }
  int Ny() const { return yend-ystart; }
  int Hx() const { return (xend-xstart)/2; }
  int Hy() const { return (yend-ystart)/2; }

  //! return true if point is inside the rectangle
  bool IsInside(const Point& Pt) const
  {
    return (Pt.x > xstart) && (Pt.x < xend) 
      &&   (Pt.y > ystart) && (Pt.y < yend);
  }
  
  bool operator==(const Window& other) const {
    return ( xstart==other.xstart && ystart==other.ystart && xend==other.xend && yend==other.yend );
  }

};

ostream& operator << (ostream& stream, const Window& w);

//! a double precision image type.
class DImage  
{

  
protected:
  int nx,ny;
  DPixel *data; 
  DPixel *data00; // address of DImage(0,0), by default same as data 
  int minindex, maxindex;
public:
  DImage(const int Nx, const int Ny);
  DImage() : data(NULL) {Allocate(0,0);}
  void Allocate(const int Nx, const int Ny, int Init=1);
  ~DImage() { if (data) delete [] data;}

  //! value on (i,j)
  //  inline DPixel& operator ()(const int i, const int j); 
  inline DPixel& operator ()(const int i, const int j) const
  {
#if DIMAGE_BOUNDS_CHECK
    if (i+j*nx<minindex || i+j*nx>maxindex) 
      {
	cerr << "CATASTROPHE : DImage out of bounds " 
	     << i << " " << j <<  " " << nx << " " << ny << " " << minindex << " " << maxindex << endl;
	abort();
      }
#endif
    return data00[i+j*nx]; /* fortran/FITS convention */
  }
  
  //! acces to pointer on first pixel of array
  DPixel* begin() const {return data;}

  DPixel* end() const {return data+nx*ny;}

  int Nx() const {return nx;}
  int Ny() const {return ny;}
  //! set all pixels to zero
  void Zero();
  //! returns the minimum pixel value 
  DPixel MinValue() const;
  //! returns the maximum pixel value 
  DPixel MaxValue() const;
  //! returns both min and max in a single image traversal. 
  void MinMaxValue(DPixel &Min, DPixel &Max) const;
  DImage(const DImage &Other);
  void dump();

  //! sum of all pixels
  DPixel sum() const;

  //! sum of pixel squares
  DPixel sum2() const;

  //! normalize the full DImage
  void Normalize();


#ifndef SWIG
  //! handy operators
  DImage& operator = (const DImage&);
  DImage& operator = (const Image&);
  DImage& operator += (const DImage &Right);
  DImage& operator -= (const DImage &Right);
  DImage& operator *= (const double &Right);
  DImage& operator *= (const DImage &Right);
  DImage& operator /= (const double &Right);
  DImage& operator /= (const DImage &Right);
  DImage& operator += (const double &Right);
  DImage& operator -= (const double &Right);
  DImage& operator = (const double &Right);
#endif

  //! write and read the DImage as a FITS array
  void readFits(const string &FitsName);
  void writeFits(const string &FitsName) const;
  //! constructor calls read routine
  DImage(const string &FitsName);
  void readFromImage(const string& FitsFileName, const Window &Rect,DPixel value_when_outside_fits=0);
  void writeInImage (const string& FitsFileName, const Window &Rect) const;
};

inline DImage operator + (const DImage &a, const DImage &b)
{
  DImage c=a;
  c += b;
  return c;
}

inline DImage operator - (const DImage &a, const DImage &b)
{
  DImage c=a;
  c -= b;
  return c;
}

inline DImage operator * (const DImage &a, const double &s)
{
  DImage c=a;
  c *= s;
  return c;
}

inline DImage operator / (const DImage &a, const double &s)
{
  DImage c=a;
  c /= s;
  return c;
}

//! An odd size DImage addressed with (0,0) at center 
//! allows quick computation of convolution like operations.
class Kernel : public DImage {

protected :
int hSizeX, hSizeY;
public:
  //! construtor for a square kernel of half size Hsize 
  Kernel(const int HSize);
  //! construtor for a rectangle kernel of half sizes HsizeX and HsizeY
  Kernel(const int HSizeX, const int HSizeY);
  Kernel(const DImage &Dim);
  Kernel(const string &FitsName);
  //! builds a bigger kernel than a given one, and sets extra pixels in the surrounding band to zero
  Kernel(const Kernel& K, int BandX, int BandY);
  Kernel(const Kernel& Other);
  Kernel();


  Kernel& operator = (const Kernel &Right); 
  void readFits(const string &FitsName);
  void readFromImage(const string& FitsFileName, const Window &Rect,DPixel value_when_outside_fits);
  bool IsDirac() const  {return (sum() == *data00);}
  //! half the size in x direction
  int HSizeX() const { return hSizeX;}
  //! half the size in y direction
  int HSizeY() const { return hSizeY;}
  //! prints all kernel values
  void dump() const ;
  //! computes mean of the kernel on (x,y)
  void bias(double &x, double &y) const;
  //! prints bias and moments of the kernel
  void dump_info(ostream &stream = cout) const;
  //! computes x and y first order moments of the kernel (\sigma_x, \sigma_y, \rho_{xy})
  void moments(double &vx, double &vy, double &vxy) const;

  bool IsEmpty() const {return (hSizeX+hSizeY == 0);}
  
  void Allocate(const int Nx, const int Ny, int Init=1);
  void allocate(const int HSizeX, const int HSizeY);
   

  friend ostream& operator << (ostream& stream, const Kernel& k);
};

//! flip the image along both axes. Does not work in place (i.e Out shout not be In).
void Mirror(const DImage &In, DImage &Out);
//! 
void Convolve(DImage& Result, const DImage& Source, const DImage &Kern);
//!
void ConvolveKernels(Kernel &Result, const Kernel &Psf, const Kernel &Kern);

#endif // DIMAGE__H
