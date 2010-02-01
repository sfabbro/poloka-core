// This may look like C code, but it is really -*- C++ -*-


#ifndef IMAGE__H
#define IMAGE__H
#include <iostream>
#include <algorithm> /* for min .. max */
//#include "persistence.h"

using namespace std;

#define BADPIX -999
#define BIGVALUE  1e+30
class Frame; /* for prototyping */



/*! \file 
    \brief basic utilities for image algebra.
*/


typedef float Pixel;

//! Class for the basic manipulation of images. 
/*!
     The internal representation uses
     (32 bits) float numbers, designated later as the Pixel type. 
*/
class Image {
  friend class FitsImageArray;
  
 public:
  //! constructor reserves the needed space to store the pixels.
      Image(const int Nx, const int Ny);

  //! reserves an empty image. 
      Image();
  //! copy constructor  
      Image(const Image&);
      virtual ~Image();

  //#define CHECK_BOUNDS

  //! access to the data (RW mode). The first pixel is indexed (0,0). 
      inline Pixel& operator()(const int i, const int j) 
	{
#ifdef CHECK_BOUNDS
	    if (i<0 || i>=nx || j<0 || j>=ny)
	      {
		cerr << "Out of range i,j,nx,ny = "<< i <<","<< j <<","<< nx <<","<< ny <<endl;
		abort();
	      }
#endif
	    return data[i+j*nx];}

  /* The following routine should be Pixel ... () const.
     It is the ImageWindow class that enforces this strange 
     signature. We should get rid of ImageWindow and PixelIterator
     and we could then have a reasonnable signature of this access 
     routine.
  */

      inline Pixel& operator() (const int i, const int j) const 
        {
#ifdef CHECK_BOUNDS
	    if (i<0 || i>=nx || j<0 || j>=ny)
	      {
		cerr << "Out of range i,j,nx,ny = "<< i <<","<< j <<","<< nx <<","<< ny <<endl;
		abort();
	      }
#endif
            return data[i+j*nx];}

      inline Pixel& value(const int i, const int j)
        { return data[i+j*nx];}

  //! returns the minimum pixel value 
      Pixel MinValue() const;
  //! returns the maximum pixel value 
      Pixel MaxValue() const;
  //! returns both min and max in a single image traversal. 
      void MinMaxValue(Pixel *Min, Pixel *Max) const;

  void ClippedMeanSigmaValue(double & Mean, double & Sigma, 
			     Image *pmask=NULL) const ;

  //! computes the median in a region of an image
  Pixel MedianInFrame(const Frame &Region, Pixel &Sigma) const;

  //! enforce the min and max value of the image.
      void EnforceMinMax(Pixel min, Pixel max) const;

      int minx() const {return 0;};
      int maxx() const {return nx-1;};
      int miny() const {return 0;};
      int maxy() const {return ny-1;};

  //! computes the mean and sigma of an image using 5 loops over ALL image. 
      void MeanSigmaValue(Pixel *Mean, Pixel *Sigma) const;
  //! computes the mean and sigma of an image using the median (from 10000 random pixels) and loop 3 times.
      void SkyLevel(Pixel *Mean, Pixel *Sigma) const;
      void SkyLevel(const Frame &Frame, Pixel *Mean, Pixel *Sigma) const;

      void dump() const;
  //! returns x size of the image. 
      int Nx() const { return nx;}
  //! returns y size of the image. 
      int Ny() const { return ny;}

#ifndef SWIG
      Image & operator =(const Image& Right);
#endif

  //! returns number of pixels 
      int NPix() const { return nx*ny;}

#ifndef SWIG
    //! 
      Image operator +(const Image& Right) const;
    //!
      Image operator -(const Image& Right) const;
    //!
      Image operator *(const Image& Right) const;
    //!
      Image operator /(const Image& Right) const;

    //!
      void operator +=(const Image& Right) const;
    //!
      void operator -=(const Image& Right) const;
    //!
      void operator *=(const Image& Right) const;
    //!
      void operator /=(const Image& Right) const;

    //!
      Image operator +(const double Right) const;
    //!
      Image operator -(const double Right) const;
    //!
      Image operator *(const double Right) const;
    //!
      Image operator /(const double Right) const;

    //!
      friend Image operator +(const double Left, const Image &Right);
    //!
      friend Image operator -(const double Left, const Image &Right);
    //!
      friend Image operator *(const double Left, const Image &Right);
    //!
      friend Image operator /(const double Left, const Image &Right);

    //!
      void operator =(const double Right);
    //!
      void operator +=(const double Right);
    //!
      void operator -=(const double Right);
    //!
      void operator *=(const double Right);
    //!
      void operator /=(const double Right);
#endif


  //! multiply by the square of the image Right.
  void MultiplyBySquare(const Image& Right);

  //! pass the image trough an heaviside function that remove negative pixels 
      void Heavyside();
  
  //! extract a subimage. 
      Image Subimage(const int x, const int y, const int width, const int height) const;

  //! extract a subimage. 
      Image Subimage(const Frame &frame) const;

  //! multiply a subimage by a double. 
      void SubimageMultiply(const int x, const int y, const int width, const int height, double factor);

  //! multiply a subimage by a double. 
      void SubimageMultiply(const Frame &frame, double factor);

  //! keep the pixels inside the mask and put the other at 0.0 
      Image Mask(const int x_Beg, const int y_Beg, const int x_End, const int y_End) const;


  //! put the pixels outside  mask to MaskValue.
      void Masking(const Frame &frame, const Pixel &MaskValue=1);

  //! put the pixels outside a disk to 1
      void DiskMaskIt(const double &xc, const double &yc, const double &radius);

  //! put the pixels outside a disk to 0
      void  Truncate(const double &xc, const double &yc, const double &radius);


  //!
      void Masking(const int x_Beg, const int y_Beg, 
		   const int x_End, const int y_End, 
		   const Pixel &MaskValue=1);

  //! put pixels above threshold to 1 and to 0 otherwise 
      void Simplify(double threshold, int above_val=1, int under_val=0);

  //! sum  pixels values
  double SumPixels(); 
  //! sum  pixels on a square centered on (x,y), side = 2*demi_size+1.
  //! not optimized but bounds are checked.
  //! used for debugging in simulation.cc
  double SumPixels(double x, double y, double demi_size) const ;


  //! sum of the squared pixel values.
  double SumSquaredPixels() const;


  void MedianFilter(const int HalfWidth);

  //! returns the pointer to the first pixel.
      Pixel* begin() { return data;}
      Pixel* begin() const { return data;}
  //! returns the pointer to the next to last pixel, as usual for containers.
      Pixel* end() { return data+nx*ny;} 
      Pixel* end() const { return data+nx*ny;} 
  //!
  bool SameSize(const Image &Other) const 
  { return (nx == Other.nx && ny == Other.ny);};
  
  //ostream& friend operator << (ostream& stream,  const Image& I)
  
  //!Builds a cosmic map
  int LaplacianFilter( const double &Sigma, const double &Mean, 
		       const double &seeing, Image &CosmicImage);
  
  //! returns the final cosmic map(after iterations) 
  //of the image using laplacian filter
  void Cosmics(const double &Sigma, const double &Mean,
	       const double &seeing, Image &CosmicImage);


  void ApplyFun(double (&F)(double))
  {
    Pixel *pend = end();
    for (Pixel *p = begin(); p < pend; ++p) *p = F(*p);
  }

  void allocate(int Nx, int Ny, int Init=1 );

 protected :


      Pixel *data;
      int nx,ny;
      Pixel *get_elem_ref(int i, int j) {return data+i+nx*j;}

 private :
      friend void image_copy(Image *To, const Image *From);


};


//! Returns the empirical covariance of 2 images.
double EmpiricCov(const Image &Im1,const Image &Im2);

//! return a clipped average of I*sqrt(W)
double ImageAndWeightError(const Image &I, const Image &W,
			 const double MinWeight = 0, double *AverageShift = NULL);  


#endif /* IMAGE__H */
