#include "pixelblock.h"


#include <iostream>


#define IMAGE_OPERATION_IMAGE(OPE)\
  if (!this->SameSize(Right))					\
    {								\
      std::cout << " PixelBlock Mismatch types in algebric operator ";	\
      abort();							 \
    }\
  PixelBlock res((const IntFrame&) (*this));		\
  const PixelType *a = this->Data();			\
  const PixelType *b = Right.Data();				\
  PixelType *r = res.Data();				\
  int size = Ntot();				\
  for (int i=size ; i>0 ; i--)			\
    {						\
      *r = *a OPE *b;				\
      ++a;					\
      ++b;					\
      ++r;					\
    }						\
  return res;					


#define IMAGE_OPERATION_EQUAL_IMAGE(OPE)\
  if (!this->SameSize(Right))					\
    {								\
      std::cout << " PixelBlock Mismatch types in algebric operator ";	\
      abort();							 \
    }\
  PixelType *a = this->Data();			\
  const PixelType *b = Right.Data();				\
  int size = Ntot();				\
  for (int i=size ; i>0 ; i--)			\
    {						\
      *a OPE *b;				\
      ++a;					\
      ++b;					\
    }						\


PixelBlock PixelBlock::operator - (const PixelBlock &Right) const
{
IMAGE_OPERATION_IMAGE(-)
}


void PixelBlock::operator -= (const PixelBlock &Right)
{
IMAGE_OPERATION_EQUAL_IMAGE(-=)
}





#include "fitsimage.h"
bool PixelBlock::ReadFromFits(const std::string &FileName,
			      const int XRef, const int YRef)
{
  FitsHeader head(FileName);
  int ntot = Ntot();
  IntFrame imageFrame(head);
  IntFrame stamp(this->Shift(XRef,YRef));
  IntFrame intersection = stamp*imageFrame;
  bool rc;
  if (intersection.Ntot() != ntot)
    {
      //Pixel is the pixel type returned by FitsHeader::read_image
      IntFrame covered(intersection.Shift(-XRef, -YRef));
      Array2D<Pixel> insidePix(covered);
      Pixel *pixData = insidePix.Data();
      rc = head.read_image(intersection.xmin, intersection.ymin, 
			   intersection.xmax , intersection.ymax, 
			   pixData);
      cout << "PixelBlock::ReadFromFits : for file " << FileName << "," << endl
	   << "the requested area is not contained in the actual image. padding with 0's " << endl;
      this->SetVal(0.);
      if (rc == 0)
	PIXEL_LOOP(insidePix, i, j)
	{
	  (*this)(i,j) = insidePix(i,j);
	}
    }
  else // the code above would work, but there is slightly faster way:
    {
      Array2D<Pixel> pixels((const IntFrame &)*this);
      Pixel *pixData = pixels.Data();
      rc = head.read_image(stamp.xmin, stamp.ymin, 
			   stamp.xmax, stamp.ymax, pixData);
      PixelType *pw = Data();
      if (rc==0) 
	for (int i = ntot; i; i--) 
	  {(*pw) = (*pixData); pw++; pixData++;}
    }
  return (rc==0);
   
}

bool PixelBlock::WriteFits(const std::string &FileName) const
{
  if (Nx()<=0 || Ny() <=0)
    {
      cout << " cannot write file " << FileName << " : nx, ny " << Nx() << ' ' << Ny() << endl;
      return false;
    }
  FitsImage im(FileName, Nx(), Ny());
  im.SetWriteAsFloat();
  int ntot = Ntot();
  Pixel *pw = im.begin();
  const PixelType *pixData = Data();
  for (int i = ntot; i; i--) (*pw++) = (*pixData++);
  return true;
}



#include <cmath>

void GaussianFill(PixelBlock &K, const double SigX, const double SigY, const double dx, const double dy)
{
  double ax = -0.5/(SigX*SigX);
  double ay = -0.5/(SigY*SigY);
  double sum = 0;
  double x1, y1;
  for (int y = K.ymin; y < K.ymax; ++y)
    for (int x = K.xmin; x < K.xmax; ++x)
      {
	x1 = x-dx ;
	y1 = y-dy ;
	double val = exp(x1*x1*ax + y1*y1*ay);
	K(x,y) = val;
	sum += val;
      }
  K *= (1./sum);
}

double ScalProd(const PixelBlock &B1, const PixelBlock &B2)
{
  if (!B1.SameSize(B2))					
    {								
      std::cout << " PixelBlock Mismatch types in ScalProd(PixelBlock,PixelBlock) ";	
      abort();
    }
  const PixelType *a = B1.Data();
  const PixelType *b = B2.Data();
  double sum = 0;
  int size = B1.Ntot();
  for (int i=size ; i ; i--)
    { 
      sum += (*a) * (*b);
      //if (i==1 or i==size) {cout << "sum *a *b " << sum << " " << *a << " " << *b << endl;}

      ++a;
      ++b;
    }

  // if ( isnan(sum) ){   B1.Write("B1.fits") ;   B2.Write("B2.fits") ; }

  return sum;
}

void PixelBlock:: Moments(double &Xmean, double &Ymean, 
			  double &VarXX , double &VarYY, double &VarXY) const
{
  double sum = 0;
  double sumx = 0;
  double sumy = 0;
  double sumxx = 0;
  double sumyy = 0;
  double sumxy = 0;
  double val;
  PIXEL_LOOP((*this),x,y)
    {
      val = (*this)(x,y);
      sum += val;
      sumx += x*val;
      sumy += y*val;
      sumxx += x*x*val;
      sumyy += y*y*val;
      sumxy += x*y*val;
    }
  if (sum == 0)
    {
      Xmean = 0;
      Ymean = 0;
      VarXX = -1;
      VarYY = -1;
      VarXY = 0;
    }
  else
    {
      Xmean = sumx/sum;
      Ymean = sumy/sum;;
      VarXX = sumxx/sum-Xmean*Xmean;
      VarYY = sumyy/sum-Ymean*Ymean;
      VarXY = sumxy/sum - Xmean*Ymean;
    }
}
