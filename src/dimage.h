// This may look like C code, but it is really -*- C++ -*-
#ifndef DIMAGE__H
#define DIMAGE__H

#include <list>
#include <string>

#include "basestar.h"
//#include "image.h"

#define DIMAGE_BOUNDS_CHECK 1

typedef double DPixel;
class Image;

//! yet another rectangle/frame/window/subimage class. used in Stamp and Vignet
struct Window { 
  Window() : xstart(0), ystart(0), xend(0), yend(0) {}
  Window(const int Xstart, const int Ystart, 
	 const int Xend, const int Yend) : xstart(Xstart), ystart(Ystart), xend(Xend), yend(Yend) {}
  int xstart, ystart, xend, yend; 
  int Nx() const { return xend-xstart; }
  int Ny() const { return yend-ystart; }
};

ostream& operator << (ostream& stream, const Window& w);

//! a double precision image type.
class DImage  
{

private :
  int nx,ny;
  DPixel *data;
protected: 
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
      }
#endif
    return data00[i+j*nx]; /* fortran/FITS convention */
  }
  
  //! acces to pointer on first pixel of array
  DPixel* begin() const {return data;}
  int Nx() const {return nx;}
  int Ny() const {return ny;}
  //! set all pixels to zero
  void Zero() { memset(data,0,sizeof(DPixel)*nx*ny);}
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
  void readFromImage(const string& FitsFileName, const Window &Rect);
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


//! an odd size DImage extracted from an Image, centered on (xc,yc)
class Stamp {  
private:
  int hsize;
public :
  //! center in the source Image
  int xc,yc ;  
  Window rect;
  //! pixels (in double precision until one changes DPixel definition )
  DImage source;  
  //! will be filled and used to discard outliers from the fit 
  double chi2; 
  int HSize () const {return hsize;}  
  //! we should not assume that there is a BaseStar corresponding to this stamp, but if any, put its pointer here
  BaseStarRef star; 
  Stamp(const double Xc, const double Yc, const Image& image, int HStampSize, BaseStar *Star);

  Stamp() :hsize(0) {};

  // need a virtual routine for Root
  virtual ~Stamp() {};
};

//! Nothing but a list of Stamps
class StampList : public list<Stamp>
{
public :
  StampList(const Image &image, const BaseStarList &starList , 
	    const int hStampSize, const int MaxStamps);

  // need a virtual routine for ROOT
  virtual ~StampList() {}
};

typedef list<Stamp>::iterator StampIterator;
typedef list<Stamp>::const_iterator StampCIterator;

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
  void readFromImage(const string& FitsFileName, const Window &Rect);
  bool IsDirac() const  {return (sum() == *data00);}
  //! half the size in x direction
  int HSizeX() const { return hSizeX;}
  //! half the size in y direction
  int HSizeY() const { return hSizeY;}
  //! zero out outside a given radius
  void KeepCircleOnly(const double &radius);
  //! fills the kernel array with a 2d normalized gaussian (xc and yc are in kernel coordinates)
  void FillWithGaussian(const double &xc, const double &yc, 
			const double &sigmax, const double &sigmay, 
			const double &rho=0);
  //! prints all kernel values
  void dump() const ;
  //! computes mean of the kernel on (x,y)
  void bias(double &x, double &y) const;
  //! prints bias and moments of the kernel
  void dump_info(ostream &stream = cout) const;
  //! computes x and y first order moments of the kernel (\sigma_x, \sigma_y, \rho_{xy})
  void moments(double &vx, double &vy, double &vxy) const;
  bool IsEmpty() const {return (hSizeX+hSizeY == 0);}
  void MaxPixel(double &xmax, double &ymax) const;
  void MinPixel(double &xmin, double &ymin) const;
  
  void Allocate(const int Nx, const int Ny, int Init=1);
   

  friend ostream& operator << (ostream& stream, const Kernel& k);
};


void Convolve(DImage& Result, const DImage& Source, const Kernel &Kern);
void ConvolveKernels(Kernel &Result, const Kernel &Psf, const Kernel &Kern);

#endif // DIMAGE__H
