#include <iostream>
#include <iomanip>
#include <math.h>

#include "dimage.h"
#include "frame.h"
#include "fileutils.h"
#include "image.h"

#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif

//***********************  DImage ***************************
DImage::DImage(const int Nx, const int Ny) : data(NULL)
{
  Allocate(Nx,Ny);
}


void DImage::Allocate(const int Nx, const int Ny, int Init)
{
  nx = Nx; ny = Ny;
  if (data) delete [] data;
  if (nx*ny) 
    {
      data = new DPixel[nx*ny];
      if (Init) Zero();
    }
  else data = NULL;
  data00 = data;
  minindex=0; maxindex = max(nx*ny-1,0);
}

DImage::DImage(const DImage& Other) : data(0)
{
  Allocate(Other.nx, Other.ny, 0);
  if (Other.data) memcpy(data, Other.data, sizeof(DPixel)*nx*ny);
}

DImage& DImage::operator = (const DImage& Right) 
{
  Allocate(Right.nx,Right.ny,0);
  if (Right.data) memcpy(data, Right.data, sizeof(DPixel)*nx*ny);
  return *this;
}
DImage& DImage::operator = (const Image& Right) 
{
  Allocate(Right.Nx(),Right.Ny(),0);
  DPixel *dp = begin();
  Pixel *p = Right.begin();
  for (int i=ny*nx; i ; --i) {*p = *dp ; ++p; ++dp;}
  return *this;
}

DPixel DImage::MinValue() const
{
  int i;
  DPixel *p;
  int npix = nx*ny;
  DPixel minval = 1e100;
  for (p=data, i=0; i<npix; ++i, ++p) if (*p < minval) minval = *p;
  return minval;
}

DPixel DImage::MaxValue() const
{
  int i;
  DPixel *p;
  int npix = nx*ny;
  DPixel maxval = -1e100;
  for (p=data, i=0; i<npix; ++i, ++p) if (*p > maxval) maxval = *p;
  return maxval;
}

void DImage::MinMaxValue(DPixel &Min, DPixel &Max) const
{
  register int i;
  register DPixel *p;
  int npix = nx*ny;
  Max = -1e100;
  Min = 1e100;
  for (p=data, i=0; i<npix; ++i, ++p)
    {
      if (*p > Max) Max = *p;
      if (*p < Min) Min = *p;
    }
}

void DImage::dump()
{
  for (int i=0; i<nx; ++i)
    {
    for (int j=0; j<ny; ++j) cout << (*this)(i,j) << ' ' ;
    cout << " endl " << endl;
    }
}

#ifdef STORAGE
DPixel& DImage::operator ()(const int i, const int j)
{ 
#ifdef DIMAGE_BOUNDS_CHECK
  if (i+j*nx<minindex || i+j*nx>maxindex)
    cerr << "CATASTROPHE : DImage out of bounds " << i << " " << j << endl;
#endif
  return data00[i+j*nx]; /* fortran/FITS convention */
}
#endif

DPixel DImage::sum() const
{
  double sum =0;
  DPixel *p = data;
  for (int i=ny*nx; i ; --i) {sum += *p ; ++p;}
  return sum;
}
void DImage::Normalize()
{
  double norm = sum();
  DPixel *p = data;
  for (int i=ny*nx; i ; --i) {*p /= norm ; ++p;}
}

DImage& DImage::operator *= (const double &Right)
{
  DPixel *p = data;
  for (int i=ny*nx; i ; --i) {*p *= Right ; ++p;}
  return *this;
}

DImage& DImage::operator += (const double &Right)
{
  DPixel *p = data;
  for (int i=ny*nx; i ; --i) {*p += Right ; ++p;}
  return *this;
}

DImage& DImage::operator += (const DImage &Right)
{
  DPixel *p = data;
  DPixel *pright = Right.data;
  for (int i=ny*nx; i ; --i) {*p += *pright; ++p; ++pright;}
  return *this;
}

DImage& DImage::operator -= (const double &Right)
{
  DPixel *p = data;
  for (int i=ny*nx; i ; --i) {*p -= Right ; ++p;}
  return *this;
}


DImage& DImage::operator -= (const DImage &Right)
{
  DPixel *p = data;
  DPixel *pright = Right.data;
  for (int i=ny*nx; i ; --i) {*p -= *pright; ++p; ++pright;}
  return *this;
}

DImage& DImage::operator = (const double &Right)
{
  DPixel *p = data;
  for (int i=ny*nx; i ; --i) {*p = Right; ++p;}
  return *this;
}

DImage& DImage::operator *= (const DImage &Right)
{
  DPixel *p = data;
  DPixel *pright = Right.data;
  for (int i=ny*nx; i ; --i) {*p *= *pright; ++p; ++pright;}
  return *this;
}

DImage& DImage::operator /= (const DImage &Right)
{
  DPixel *p = data;
  DPixel *pright = Right.data;
  for (int i=ny*nx; i ; --i) {*p /= *pright; ++p; ++pright;}
  return *this;
}


#include "fitsimage.h"
void DImage::writeFits(const string &Name) const
{
  Image fdata(nx,ny);
  DPixel *dp = begin();
  Pixel *p = fdata.begin();
  for (int i=ny*nx; i ; --i) {*p = *dp ; ++p; ++dp;}
  FitsImage fits(Name, fdata); 
  fits.SetWriteAsFloat();
}

DImage::DImage(const string &FitsName) : data(NULL)
{
  readFits(FitsName);
}

void DImage::readFits(const string &FitsName)
{
  if (!FileExists(FitsName))
    {
      Allocate(0,0);
      return;
    }
  FitsImage im(FitsName);
  Allocate(im.Nx(),im.Ny());
  DPixel *dp = begin();
  Pixel *p = im.begin();
  for (int i=ny*nx; i ; --i) {*dp = *p ; ++p; ++dp;}
}


//***********************  Stamp ***************************

Stamp::Stamp(const double Xc, const double Yc, const Image& image, int HStampSize, BaseStar *Star ) 
  : hsize(HStampSize), xc(int(Xc)), yc(int(Yc)), source(0,0), star(Star)
{
  // cout << " debug xc yc " << xc << ' ' << yc << endl;
  int width = 2*HStampSize+1;
  int height = 2*HStampSize+1;
  source.Allocate(width, height);
  int xoff = xc - HStampSize; 
  int yoff = yc - HStampSize;
  xstart = max(0,xc - HStampSize);
  ystart = max(0,yc - HStampSize);
  int xend = min(xstart+width,image.Nx());
  int yend = min(ystart+height, image.Ny());
  if (width*height != (xstart-xend)*(ystart-yend)) 
    cerr << " we miss pixels for (" << xc << "," << yc << ")" << endl;
  for (int j=ystart; j <yend; ++j)
    {
      for (int i= xstart; i< xend; i++)
	{
	  source(i-xoff,j-yoff) = image(i,j);
	}
    }
  //  cout << " debug extraction " <<  image(xc,yc) << ' ' << source(HStampSize,HStampSize) << " xstart y " << xstart << ' ' << ystart << endl;
}

StampList::StampList(const Image &image, const BaseStarList &starList, const int hStampSize, const int MaxStamps)
{
  int count = 0;
  Frame imageFrame(Point(0,0),Point(image.Nx(), image.Ny()));
  for (BaseStarCIterator si = starList.begin(); si != starList.end() && count < MaxStamps; ++si)
    {
      const BaseStar *s = *si;
      if (imageFrame.MinDistToEdges(*s) < hStampSize+2) continue;
      Stamp stamp(s->x, s->y, image, hStampSize, (BaseStar *) s);
      push_back(stamp); 
      ++count;
    }
}


//*********************** Kernel ***************************

Kernel::Kernel(const int HSizeX, const int HSizeY) :  
  DImage(2*HSizeX+1, 2*HSizeY+1), hSizeX(HSizeX), hSizeY(HSizeY)
{
  data00 = &(*this)(hSizeX,hSizeY);
  minindex = begin()-data00; maxindex = minindex+Nx()*Ny()-1;
};

Kernel::Kernel(const int HSize) :
  DImage(2*HSize+1, 2*HSize+1), hSizeX(HSize), hSizeY(HSize)
{
  data00 = &(*this)(hSizeX,hSizeY); 
  minindex = begin()-data00; maxindex = minindex+Nx()*Ny()-1;
};

Kernel::Kernel() : DImage(), hSizeX(0), hSizeY(0) {}

Kernel::Kernel(const Kernel& Other) : DImage(Other) 
{ 
  hSizeX = Other.hSizeX; 
  hSizeY = Other.hSizeY; 
  data00 = &(*this)(hSizeX,hSizeY);
  minindex = begin()-data00; 
  maxindex = minindex + Nx()*Ny() - 1;
}

Kernel& Kernel::operator = (const Kernel &Right)
{ 
  *((DImage*) this)  =  Right; 
  hSizeX = Right.hSizeX; 
  hSizeY = Right.hSizeY; 
  data00 = &(*this)(hSizeX,hSizeY);
  minindex = begin()-data00; 
  maxindex = minindex + Nx()*Ny()-1;
  return *this;
}

Kernel::Kernel(const DImage &Dim) : DImage(Dim)
{  
  hSizeX = (Nx()-1)/2;
  hSizeY = (Ny()-1)/2;
  data00 = &(*this)(hSizeX,hSizeY);
  minindex = begin()-data00; 
  maxindex = minindex + Nx()*Ny()-1;
}

Kernel::Kernel(const string &FitsName) : DImage(FitsName)
{  
  hSizeX = (Nx()-1)/2;
  hSizeY = (Ny()-1)/2;
  data00 = &(*this)(hSizeX,hSizeY);
  minindex = begin()-data00; 
  maxindex = minindex + Nx()*Ny()-1;
}

void Kernel::readFits(const string &FitsName)
{
  FitsImage im(FitsName);
  if ((im.Nx() == 0) && (im.Ny() == 0) ) return;
  DImage::readFits(FitsName);
  hSizeX = (Nx()-1)/2;
  hSizeY = (Ny()-1)/2;
  data00 = &(*this)(hSizeX,hSizeY);
  minindex = begin()-data00; 
  maxindex = minindex + Nx()*Ny()-1;
}

void Kernel::dump() const
{
  if (maxindex <= minindex) return;
  for (int j=-hSizeY; j<= hSizeY ; ++j) 
    {
      for (int i=-hSizeX; i<=hSizeX;++i) printf("%7.4f ",(*this)(i,j));
      printf("\n");
    }
}

void Kernel::bias(double &x, double &y) const
{
  if (maxindex <= minindex) return;
  double sumx=0;
  double sumy=0;
  double sum =0;
  for (int j=-hSizeY; j <= hSizeY ; ++j)
    for (int i=-hSizeX; i <= hSizeX  ;++i)
      {
	double value = (*this)(i,j);
	sumx += i*value;
	sumy += j*value;
	sum += value;
      }
  if (sum) {x = sumx/sum; y = sumy /sum; }
 else {x = -1000.; y = -1000;}
}


void Kernel::moments(double &vx, double &vy, double &vxy) const
{
  if (maxindex <= minindex) return;
  double sum=0;
  double sumx = 0;
  double sumy=0;
  double sumx2=0;
  double sumy2=0;
  double sumxy=0;
  for (int j=-hSizeY; j <= hSizeY ; ++j)
    for (int i=-hSizeX; i <= hSizeX  ;++i)
      {
       double value = (*this)(i,j);
       sum += value;
       sumx += i*value;
       sumy += j*value;
       sumx2 += i*i*value;
       sumy2 += j*j*value;
       sumxy += i*j*value;
     }
 sumx /= sum;
 sumy /= sum;
 vx = sumx2/sum - sumx*sumx;
 vy = sumy2/sum - sumy*sumy;
 vxy = sumxy/sum - sumx*sumy;
}


void Kernel::dump_info(ostream &stream) const
{
  if (maxindex <= minindex) return;
  double dx,dy, varx, vary, varxy;
  bias(dx,dy);
  moments(varx,vary,varxy);
  int oldprec = stream.precision();
  stream << setprecision(10);
  stream << " bias ( " << dx << ',' << dy << 
    ") sx,sy,rho, sum " << sqrt(varx) << ',' << sqrt(vary) << ',' << varxy/sqrt(varx*vary) << ',' << sum () << endl;
  stream <<setprecision(oldprec);
}

Kernel::Kernel(const Kernel& K, int BandX, int BandY) 
  : DImage(K.Nx()+2*BandX, K.Ny()+2*BandY)
{  
  hSizeX = (Nx()-1)/2;
  hSizeY = (Ny()-1)/2;
  data00 = &(*this)(hSizeX,hSizeY);
  minindex = begin()-data00; 
  maxindex = minindex + Nx()*Ny()-1;
  // copy only common part with K
  int hkx = K.HSizeX();
  int hky = K.HSizeY();
  int kx = (hSizeX > hkx) ? hkx : hSizeX;
  int ky = (hSizeY > hky) ? hky : hSizeY;
  for (int j=-ky; j <= ky ; ++j)
    for (int i=-kx; i <= kx  ;++i)
      (*this)(i,j) = K(i,j);

}
void Kernel::MaxPixel(double &xmax, double &ymax)
{
  double maxval = -1e30;
  for (int j=-hSizeY; j <= hSizeY ; ++j)
    for (int i=-hSizeX; i <= hSizeX  ;++i)
      {
	double value = (*this)(i,j);
	if (value > maxval) {maxval=value;xmax=i;ymax=j;}
      }
}

void Kernel::MinPixel(double &xmin, double &ymin)
{
  double minval = 1e30;
  for (int j=-hSizeY; j <= hSizeY ; ++j)
    for (int i=-hSizeX; i <= hSizeX  ;++i)
      {
	double value = (*this)(i,j);
	if (value < minval) {minval=value;xmin=i;ymin=j;}
      }
}

void Kernel::KeepCircleOnly(const double &radius)
{
  if (maxindex <= minindex) return;
  double rad2  = radius*radius;
  for (int j=-hSizeY; j <= hSizeY ; ++j)
    for (int i=-hSizeX; i <= hSizeX  ;++i)
      if ( i*i + j*j > rad2) (*this)(i,j) = 0;
}

void Kernel::FillWithGaussian(const double &xc, const double &yc, 
			      const double &sigmax, const double &sigmay, 
			      const double &rho)
{
  double beta = 1-rho*rho;
  double alphax = -0.5/(sigmax*sigmax*beta);
  double alphay = -0.5/(sigmay*sigmay*beta);
  double alphaxy = rho/(sigmax*sigmay);
  double norm = 1./(2.*M_PI*sigmax*sigmay*sqrt(beta));
  for (int j=-hSizeY; j<=hSizeY; ++j)
    for (int i=-hSizeX; i<=hSizeX; ++i)
      {
	double x = double(i-xc);
	double y = double(j-yc);
	(*this)(i,j) = exp(x*x*alphax + y*y*alphay + x*y*alphaxy)*norm;
      }
}


#define OPTIMIZED

void Convolve(DImage& Result, const DImage& Source, const Kernel &Kern)
{
 /* assumes that Result is already allocated and properly sized */
#ifdef OPTIMIZED
  int nstepx = Kern.Nx();
#endif
  int endrx = Result.Nx();
  int endry = Result.Ny();
  int ksx = Kern.HSizeX();
  int ksy = Kern.HSizeY();
 
  for (int j=0; j< endry; ++j)
    for (int i=0; i < endrx; ++i)
    {
      double sum = 0;
#ifndef OPTIMIZED  
      for (int jk = -ksy; jk <= ksy; ++jk)
	for (int ik = -ksx; ik <= ksx; ++ik)
	  sum += Kern(ik,jk)*Source(i+ksx-ik,j+ksy-jk);
#else
      DPixel *pk = Kern.begin();
      for (int jk = -ksy; jk <=ksy; ++jk)
	{
	  DPixel *ps = &Source(i+2*ksx,j+ksy-jk);
	  for (int toto = nstepx; toto; --toto) {sum += (*pk) * (*ps); ++pk; --ps;}
	}
#endif
      Result(i,j) = sum;
    }
}

void ConvolveKernels(Kernel &Result, const Kernel &Psf, const Kernel &Kern)
{
  // you can also allocate Result before hand
  // recall : Result.HSize() = Psf.HSize() - Kern.HSize(), outside of it it is undefined

  int hpx = Psf.HSizeX();
  int hpy = Psf.HSizeY();
  int hkx = Kern.HSizeX();
  int hky = Kern.HSizeY();
  int hrx = hpx-hkx;
  int hry = hpy-hky;
  if ((Result.HSizeX() != hrx) || (Result.HSizeY() != hry)) 
    { 
      cout << " Changing result size in ConvolveKernels: " << endl;
      cout << "   hx : " << Result.HSizeX() << " -> " << hrx << endl
	   << "   hy : " << Result.HSizeY() << " -> " << hry << endl;
      Result = Kernel(hrx,hry);
    }

  for (int j=-hry; j<=hry; ++j)
    for (int i=-hrx; i<=hrx; ++i)
      {
	double sum = 0;
#ifndef OPTIMIZED  
	for (int jk =-hky; jk <= hky; ++jk)
	  for (int ik =-hkx; ik <= hkx; ++ik)
	    sum += Kern(ik,jk) * Psf(i-ik,j-jk);
#else
	DPixel *pkern = Kern.begin();
	for (int jk =-hky; jk <= hky; ++jk)
	  {
	    DPixel *ppsf = &Psf(i+hkx, j-jk);
	    for (int ik = -hkx; ik <= hkx; ++ik) 
	      {sum += (*pkern) * (*ppsf); ++pkern; --ppsf;}
	  }
#endif
	Result(i,j) = sum;
      }
}


#ifdef USE_ROOT
/*
RUN_ROOTCINT
LINKDEF_CONTENT : #pragma link C++ class DImage+;
LINKDEF_CONTENT : #pragma link C++ class Stamp+;
LINKDEF_CONTENT : #pragma link default on;
LINKDEF_CONTENT : #pragma link C++ class list<Stamp>;
LINKDEF_CONTENT : #pragma link C++ class StampList+;
LINKDEF_CONTENT : #pragma link C++ class Kernel+;

*/
#include "root_dict/dimagedict.cc"

#endif /* USE_ROOT */



