#include <iostream>
#include <iomanip>
#include <string.h>
#include <algorithm>
#include <math.h>

#include "image.h"
#include "imageback.h"
#include "frame.h"
#include "vutils.h"



Image::Image(const int Nx, const int Ny) : data(NULL)
  {
  allocate(Nx, Ny);
  }

Image::Image()
{
nx=ny=0;
data = NULL;
}

Image::Image(const Image& Source) : data(NULL)
  {
   image_copy(this, &Source);
  }

void image_copy(Image *To, const Image *From)
{
To->allocate(From->nx, From->ny, 0);
memcpy(To->data, From->data, From->nx * From->ny * sizeof(Pixel));
}


Image & Image::operator =(const Image& Right)
{
image_copy(this, &Right);
return (*this);
}



/****************   Binary Operations with Images *************/

static bool same_sizes(const Image& I1, const Image& I2)
{
  if ((I1.Nx() == I2.Nx()) && (I1.Ny() == I2.Ny())) return true;
  cerr << " trying to perform an operation on images with different sizes" << endl;
  abort(); // we should in fact throw an exception.
  return false;
}



#define IMAGE_OPERATION_IMAGE(OPE)\
{\
same_sizes(*this,Right);\
Image res(this->nx, this->ny);\
Pixel *a = this->data;\
Pixel *b = Right.data;\
Pixel *r = res.data;\
int size = this->nx * this->ny;\
for (int i=size ; i>0 ; i--)\
  {\
    *r = *a OPE *b;\
    ++a;\
    ++b;\
    ++r;\
  }\
return res;\
}\

Image Image::operator +(const Image& Right) const
{
IMAGE_OPERATION_IMAGE(+);
}

Image Image::operator -(const Image& Right) const
{
IMAGE_OPERATION_IMAGE(-);
}

Image Image::operator *(const Image& Right) const
{
IMAGE_OPERATION_IMAGE(*);
}

Image Image::operator /(const Image& Right) const
{
IMAGE_OPERATION_IMAGE(/);
}

#undef IMAGE_OPERATION_IMAGE

#define IMAGE_OPERATION_EQUAL_IMAGE(OPE)\
{\
same_sizes(*this,Right);\
Pixel *a = this->data;\
Pixel *b = Right.data;\
int size = this->nx * this->ny;\
for (int i=size ; i>0 ; i--)\
  {\
    *a OPE *b;\
    ++a;\
    ++b;\
  }\
}\

void Image::operator +=(const Image& Right) const
{
IMAGE_OPERATION_EQUAL_IMAGE(+=);
}

void Image::operator -=(const Image& Right) const
{
IMAGE_OPERATION_EQUAL_IMAGE(-=);
}

void Image::operator *=(const Image& Right) const
{
IMAGE_OPERATION_EQUAL_IMAGE(*=);
}

void Image::operator /=(const Image& Right) const
{
IMAGE_OPERATION_EQUAL_IMAGE(/=);
}

#undef IMAGE_OPERATION_EQUAL__IMAGE

/*********************** Operations with an Image and a scalar ***************/

#define IMAGE_OPERATION_SCALAR(OPE)\
{\
Image res(this->nx, this->ny);\
Pixel *a = this->data;\
Pixel *r = res.data;\
int size = this->nx * this->ny;\
for (int i=size ; i>0 ; i--)\
  {\
    *r = *a OPE Right;\
    ++a;\
    ++r;\
  }\
return res;\
}\

Image Image::operator +(const double Right) const
{
IMAGE_OPERATION_SCALAR(+);
}

Image Image::operator -(const double Right) const
{
IMAGE_OPERATION_SCALAR(-);
}

Image Image::operator *(const double Right) const
{
IMAGE_OPERATION_SCALAR(*);
}

Image Image::operator /(const double Right) const
{
IMAGE_OPERATION_SCALAR(/);
}

#undef IMAGE_OPERATION_SCALAR



#define SCALAR_OPERATION_IMAGE(OPE)\
{\
Image res(Right.nx, Right.ny);\
Pixel *a = Right.data;\
Pixel *r = res.data;\
int size = Right.nx * Right.ny;\
for (int i=size ; i>0 ; i--)\
  {\
    *r = Left OPE *a;\
    ++a;\
    ++r;\
  }\
return res;\
}\

Image operator +(const double Left, const Image &Right)
{
SCALAR_OPERATION_IMAGE(+);
}

Image operator -(const double Left, const Image &Right)
{
SCALAR_OPERATION_IMAGE(-);
}

Image operator *(const double Left, const Image &Right)
{
SCALAR_OPERATION_IMAGE(*);
}

Image operator /(const double Left, const Image &Right)
{
SCALAR_OPERATION_IMAGE(/);
}
#undef SCALAR_OPERATION_IMAGE



#define IMAGE_OPERATION_EQUAL_SCALAR(OPE)\
{\
Pixel *a = this->data;\
int size = this->nx * this->ny;\
for (int i=size ; i>0 ; i--)\
  {\
    *a OPE Right;\
    a++;\
  }\
}\

void Image::operator =(const double Right)
{
IMAGE_OPERATION_EQUAL_SCALAR(=);
}

void Image::operator +=(const double Right)
{
IMAGE_OPERATION_EQUAL_SCALAR(+=);
}

void Image::operator -=(const double Right)
{
IMAGE_OPERATION_EQUAL_SCALAR(-=);
}

void Image::operator *=(const double Right)
{
IMAGE_OPERATION_EQUAL_SCALAR(*=);
}

void Image::operator /=(const double Right)
{
IMAGE_OPERATION_EQUAL_SCALAR(/=);
}

#undef IMAGE_OPERATION_EQUAL_SCALAR


void Image::MultiplyBySquare(const Image& Right)
{
  same_sizes(*this,Right);
  Pixel *a = this->data;
  Pixel *b = Right.data;
  int size = this->nx * this->ny;
  for (int i=size ; i>0 ; i--)
    {
       *a *= (*b)*(*b);
       ++a;
       ++b;
    }
}







void Image::Heavyside()
{
Pixel *a = this->data;
int size = this->nx * this->ny;
for (int i=size ; i>0 ; i--)
  {
    if (*a < 0) *a = 0.0;
    ++a;
  }
}
#include "gtransfo.h"

#ifdef STORAGE
double Image::BilinearInterpolate(const double & x, const double & y, int & xp, int & yp) const

{
  double dx,dy;
  // out of image ?? assumed to be already tested.
  
  dx = x - xp;
  dy = y - yp;

  /* on ecrabouille les valeurs precedentes de xp,yp,dx,dy*/
  xp = int(floor(x));
  yp = int(floor(y));
  
  dx = x - xp;
  dy = y - yp;
  
  if (xp<0){xp = 0; dx -= 1.0;}
      else if (xp>=(nx-1)) {xp = nx - 2; if (nx<2) xp = 0; dx += 1.0;}
  if (yp<0){yp = 0; dy -= 1.0;}
  else if (yp>=(ny-1)) {yp = ny - 2; if (ny<2) yp = 0; dy += 1.0;}
  
  return (((*this)(xp,yp)*(1-dx) + (*this)(xp+1,yp)*dx)*(1-dy) +
	  ((*this)(xp,yp+1)*(1-dx) + (*this)(xp+1,yp+1)*dx)*dy);
}


static Pixel QuadraticInterpolate_Old(const Image &In, 
				  const double & x, const double & y, 
				  int & xp, int & yp)

{
  // out of image ??  assumed to be already tested
  
  double dx = x -xp;
  double dy = y -yp;

  double  z00 =    (dx-1)*(dy-1)*dx*dy;
  double  z01 = -2*(dx-1)*(dy+1)*(dy-1)*dx;
  double  z02 =    (dx-1)*(dy+1)*dx*dy;
  double  z10 = -2*(dx+1)*(dx-1)*(dy-1)*dy;
  double  z11 =  4*(dx+1)*(dx-1)*(dy+1)*(dy-1);
  double  z12 = -2*(dx+1)*(dx-1)*(dy+1)*dy;
  double  z20 =    (dx+1)*(dy-1)*dx*dy;
  double  z21 = -2*(dx+1)*(dy+1)*(dy-1)*dx;
  double  z22 =    (dx+1)*(dy+1)*dx*dy;

    double interp = 
        In(xp-1, yp-1)*z00
      + In(xp  , yp-1)*z10
      + In(xp+1, yp-1)*z20
      + In(xp-1, yp  )*z01
      + In(xp  , yp  )*z11
      + In(xp+1, yp  )*z21
      + In(xp-1, yp+1)*z02
      + In(xp  , yp+1)*z12
      + In(xp+1, yp+1)*z22;
    return interp*0.25;
}
#endif /* STORAGE */


// checked that this routine does the same as the old one
// it is also slightly faster and can accomodate variance maps
static Pixel QuadraticInterpolate(const Image &In,
				  const double & x, const double & y, 
				  int & xp, int & yp, bool IsVarianceMap = false)
{

  // out of image ??  assumed to be already tested
  
  double dx = x -xp;
  double dy = y -yp;
  
  double xf0 = (dx-1)*dx;
  double xf1 = (dx-1)*(dx+1);
  double xf2 = (dx+1)*dx;
  double yf0 = (dy-1)*dy;
  double yf1 = (dy-1)*(dy+1);
  double yf2 = (dy+1)*dy;
  double minus2 = -2.;
  double one_fourth = 0.25;

  if (IsVarianceMap) // coefficients for variances are to be squared !
    {
      xf0 *= xf0;
      xf1 *= xf1;
      xf2 *= xf2;
      yf0 *= yf0;
      yf1 *= yf1;
      yf2 *= yf2;
      minus2*= minus2;
      one_fourth *= one_fourth;
    }

       double interp = 
	 (    In(xp-1, yp-1)  *xf0
	   +  In(xp  , yp-1)*xf1*minus2
	   +  In(xp+1, yp-1)*xf2) * yf0 
	 +   (In(xp-1, yp  )*xf0
	    + In(xp  , yp  )*xf1*minus2
	    + In(xp+1, yp  )*xf2) * minus2 * yf1
	 +   (In(xp-1, yp+1)*xf0
	    + In(xp  , yp+1)*xf1*minus2
	    + In(xp+1, yp+1)*xf2) * yf2;
    return interp*one_fourth;
}

Pixel Image::Interpolate(const double x, const double y, const int level,
			 const bool IsVarianceMap) const
{
  int xp,yp;

  /* closest integers, not integer part */
  xp = int(floor(x+0.5));
  yp = int(floor(y+0.5));

  if (xp<0 || yp<0 || xp > nx || yp > ny) return 0.0;


  /* if on the sides, run the bilinear interpolation (4 pixels) which handles
     side effects */

  if ((level==2) || (xp < 1) || (yp < 1) || (xp > nx-2) || (yp > ny-2))
    {
      xp = int(floor(x));
      yp = int(floor(y));
  
      double dx = x - xp;
      double dy = y - yp;
  
      // handle side effects
      if (xp<0){xp = 0; dx -= 1.0;}
      else if (xp>=(nx-1)) {xp = nx - 2; if (nx<2) xp = 0; dx += 1.0;}
      if (yp<0){yp = 0; dy -= 1.0;}
      else if (yp>=(ny-1)) {yp = ny - 2; if (ny<2) yp = 0; dy += 1.0;}
      if (!IsVarianceMap)
	{
	  return (((*this)(xp,yp)*(1-dx) + (*this)(xp+1,yp)*dx)*(1-dy) +
		  ((*this)(xp,yp+1)*(1-dx) + (*this)(xp+1,yp+1)*dx)*dy);
	}
      else // transforming variances,
	{ // the same as above with squared coefficients:
	  return (
	   ((*this)(xp,yp)*(1-dx)*(1-dx) + 
	    (*this)(xp+1,yp)*dx*dx)
	   *(1-dy)*(1-dy) 
	   +((*this)(xp,yp+1)*(1-dx)*(1-dx) + 
	     (*this)(xp+1,yp+1)*dx*dx)
	   *dy*dy);
	}
    }
  else // Quadratic interpolate
    {
      return QuadraticInterpolate(*this, x, y, xp, yp, IsVarianceMap);
    }
}


Image Image::GtransfoImage(const Gtransfo & g, int nx, int ny, 
			   float DefaultVal, const int interpLevel,
			   const bool IsVarianceMap) const
{
  
  if (IsIdentity(&g)) return *this;
  double endXStart = double(this->Nx()-1);
  double endYStart = double(this->Ny()-1);
  
  if (nx==0) nx=this->Nx();
  if (ny==0) ny=this->Ny();
  
  Image result(nx,ny);
  
  double averageJacobian = fabs(g.Jacobian(int(nx/2),int(ny/2)));

  for (int j=0; j<ny; j++) 
    for (int i=0; i<nx; i++)
      {
	double xout,yout;
	g.apply(i,j,xout,yout);
	if (xout >= 0.0 && xout < endXStart && yout >= 0.0 && yout < endYStart)
	  result(i,j) = this->Interpolate(xout,yout,interpLevel,IsVarianceMap) 
	    * averageJacobian;

	/* there is a subbtle point here : The variation of the jacobian due 
	   to the fact that pixels may have different sizes on the sky across 
	   the field was taken into account by flatfielding. 
	   taking local Jacobian reintroduces the (differential) pixel size 
	   variation */

	else
	  result(i,j) = DefaultVal * averageJacobian;
      }
  return result;
}

Image Image::Subimage(const int x, const int y, const int width, const int height) const
{
ImageWindow w(*this,x,y,width,height);
int sizex = w.width();
int sizey = w.height();
cout << " sizex " << sizex << " sizey " << sizey << endl;
Image vignette(sizex, sizey);

PixelIterator p(w);
for (int j=0; j<sizey; j++)
  {
  for (int i= 0; i< sizex; i++)
    {
    vignette(i,j) = *p;
    ++p;
    }
  }
return vignette;
}


Image Image::Subimage(const Frame &frame) const
{
  int x = (int) frame.xMin ;
  int y = (int) frame.yMin ;
  int X = (int) frame.xMax ;
  int Y = (int) frame.yMax ;
  int width = X - x ;
  int height = Y - y ;
  return(this->Subimage(x,y,width,height));
}

void Image::SubimageMultiply(const int x, const int y, const int width, const int height, double factor)
{
ImageWindow w(*this,x,y,width,height);
int sizex = w.width();
int sizey = w.height();
 cout << " Multiply subimage of sizex " << sizex << " and sizey " << sizey << " by factor : "<< factor << endl;

PixelIterator p(w);
for (int i=0; i<sizex*sizey; i++)
  {
    *p *= factor;
    ++p;
  }
}

void Image::SubimageMultiply(const Frame &frame, double factor)
{
  int x = (int) frame.xMin ;
  int y = (int) frame.yMin ;
  int X = (int) frame.xMax ;
  int Y = (int) frame.yMax ;
  int width = X - x ;
  int height = Y - y ;
  this->SubimageMultiply(x,y,width,height,factor);
}


void Image::Masking(const Frame &frame, const Pixel &MaskValue)
{
  this->Masking(int(frame.xMin), int(frame.yMin), int(frame.xMax), int(frame.yMax), MaskValue);
}


void Image::Masking(const int x_Beg, const int y_Beg, 
		    const int x_End, const int y_End, 
		    const Pixel &MaskValue) 
{

  int ny = Ny(), nx = Nx();
  cerr << " Masking : xbg, ybeg, xend, yend " << x_Beg << " " <<  y_Beg 
       << " " << x_End << " " << y_End << endl ;
  for (int j=0; j<ny; ++j) for (int i=0; i<nx; ++i)
    {
      if (j<y_Beg || j>y_End || i<x_Beg || i>x_End)
	{
	  (*this)(i,j) = MaskValue;
	}
    }
}



Image Image::Mask(const int x_Beg, const int y_Beg, const int x_End, const int y_End) const
{
int ny = Ny(), nx = Nx();
Image masked(nx,ny);

for (int j=0; j<ny; j++) for (int i=0; i<nx; i++)
  {
    if (j<y_Beg || j>y_End || i<x_Beg || i>x_End)
      {
	masked(i,j) = 0.0;
      }
    else
      {
	masked(i,j) = (*this)(i,j);
      }
  }

return masked;
}

void Image::DiskMaskIt(const double &xc, const double &yc, const double &radius)
{
  double rad2=radius*radius;
  for (int j=0; j<ny; j++) for (int i=0; i<nx; i++)
    {
      if ((j-yc)*(j-yc)+(i-xc)*(i-xc) < rad2)
	  (*this)(i,j) = 0.0;
    }
}

void Image::Truncate(const double &xc, const double &yc, const double &radius)
{
  double rad2=radius*radius;
  for (int j=0; j<ny; j++) for (int i=0; i<nx; i++)
    {
      if ((j-yc)*(j-yc)+(i-xc)*(i-xc) > rad2)
	  (*this)(i,j) = 0.;
    }
}

void Image::allocate(int Nx, int Ny, int Init)
{
if (data) delete [] data;
nx = Nx; ny = Ny;
if (nx*ny) 
  {
  data = new Pixel [nx*ny];
  //  cout << "Image::allocate data = "<< data << endl ;
  if (Init) memset(data, 0, nx*ny*sizeof(Pixel));
  }
}

Image::~Image()
  {
    //  cout << "Image::~Image data = "<< data << endl ;
  delete []data;
  }

Pixel Image::MinValue() const
{ /* works for an empty image */
int i;
Pixel *p;
int npix = nx*ny;
Pixel minval = 1e30;
for (p=data, i=0; i<npix; i++, p++) if (*p < minval) minval = *p;
return minval;
}

Pixel Image::MaxValue() const
{ /* works for an empty image */
int i;
Pixel *p;
int npix = nx*ny;
Pixel maxval = -1e30;
for (p=data, i=0; i<npix; i++, p++) if (*p > maxval) maxval = *p;
return maxval;
}


double Image::SumSquaredPixels() const
{
  double S2 =0;
  Pixel *pend = end();
  for (Pixel *p=data;  p<pend; ++p)
    {
      S2 += (*p) * (*p);
    }
  return S2;
}






void Image::MinMaxValue(Pixel *Min, Pixel *Max) const
{
register int i;
register Pixel *p;
int npix = nx*ny;
register Pixel maxval = -1e30;
register Pixel minval = 1e30;
for (p=data, i=0; i<npix; i++, p++)
  {
  if (*p > maxval) maxval = *p;
  if (*p < minval) minval = *p;
  }
*Min = minval; *Max = maxval;
}

double
Image::SumPixels() 
{
  register int i;
  register Pixel *p;
  int npix = nx*ny;
  double S = 0. ;
  for (p=data, i=0; i<npix; i++, p++)
    {
      S += *p ;
    }
  return(S);
}


// This is used in debug in simulation.cc
double Image::SumPixels(double x, double y, double demi_size) const
{
  int xmin = (int) ( x - demi_size + 0.5 );
  int ymin = (int) ( y - demi_size + 0.5 );
  int xmax = (int) ( x + demi_size + 1. + 0.5 );
  int ymax = (int) ( y + demi_size + 1. + 0.5 );
  if (xmin < 0 ) xmin = 0 ;  
  if (ymin < 0 ) ymin = 0 ;
  if (xmax > Nx() ) xmax = Nx();
  if (ymax > Ny() ) ymax = Ny();
  double S = 0. ;
  for(int i = xmin ; i < xmax ; i++)    
    for(int j = ymin ; j < ymax ; j++)
      {
	S += (*this)(i,j);
      }
  return S ;
}

void 
Image::Simplify(double threshold, int above_val, int  under_val ) 
{
  register Pixel *p, *pend;
  pend = data + nx*ny;
  for (p=data; p < pend; p++)
    {
      if (*p >threshold ) 
	*p=above_val;
      else
	*p =under_val ;
  }
}

void Image::EnforceMinMax(Pixel min, Pixel max)const
{
register int i;
register Pixel *p;
int npix = nx*ny;

for (p=data, i=0; i<npix; i++, p++)
  {
  if (*p > max)  (*p) = max;
  if (*p < min)  (*p) = min;
  }
}


void Image::ClippedMeanSigmaValue(double & Mean, double & Sigma, 
				    Image *pmask) const
{

  Mean = 0 ;
  Sigma = 0 ;

  int nx = Nx();
  int ny = Ny();

  int ntot = nx*ny;
  int npix = 50000 ;
  int pas = (int) (ntot/(1.*npix));
  if( npix > ntot )
    {
      npix = ntot;
      pas = 1 ;
    }

  
  register Pixel *p;
  register Pixel *pend = data + ntot;
  register Pixel *pm;

  double average=0;
  double rms = 1e30;

  for(int loop=0; loop<10; loop++)
    {
      int count=0;
      double meanval = 0;
      double sigval = 0;
      for (p=data, pm = pmask->data ; p<pend; p += pas, pm += pas)
	{
	  if (pmask &&  *pm > 1e-10) continue; // ignore pixel if masked
	  double val = *p - average;
	  if (fabs(val)<3.*rms)
	    {
	      meanval += val;
	      sigval += val*val;
	      count++;
	    }
	}
      if (count == 0)
	{
	  cerr << "Erreur calcul median et sigma, count=0" << endl ;
	  return ;
	}
      meanval /= count;
      average += meanval;
      sigval = (sigval/(count)) - (meanval*meanval);
      if (sigval > 0) 
	{
	  sigval = sqrt(sigval);
	  rms = sigval;
	}
    }

  Mean = average ;
  Sigma = rms;
}

Pixel Image::MedianInFrame(const Frame &Region, Pixel &Sigma) const
{
  int npix = 0;
  int x0 = max(int(Region.xMin),0);
  int y0 = max(int(Region.yMin),0);
  int x1 = min(int(Region.xMax),nx);
  int y1 = min(int(Region.yMax),ny);
  Pixel *regionArray = new Pixel[int(Region.Nx()*Region.Ny())];
  for (int j=y0; j<y1; j++)
    for (int i=x0; i<x1; i++)
      {
	regionArray[npix] = (*this)(i,j);
	npix++;
      }
  Pixel median = Fmedian_sigma(regionArray,npix,Sigma);
  delete [] regionArray;

  return median;
}



void Image::MeanSigmaValue(Pixel *Mean, Pixel *Sigma) const
{
register int i;
register Pixel *p;
int npix = nx*ny;
double meanval,average;
double sigval,rms;

average = 0.0;
rms = 1e30;

int nloop = 5;
for(int loop=0; loop<nloop; loop++)
  {
    meanval = 0.0;
    sigval = 0.0;
    int count = 0;
    for (p=data, i=0; i<npix; i++, p++)
      {
        double val = *p - average;
	if (fabs(val)<3.*rms)
	  {
	    meanval += val;
	    sigval += val*val;
	    count++;
	  }
      }
    meanval /= count;
    average += meanval;
    sigval = (sigval/(count)) - (meanval*meanval);
    if (sigval > 0) 
      {
	sigval = sqrt(sigval);
	rms = sigval;
      }
    // cerr << "count : " << count << endl ; 
  }
*Mean = average; *Sigma = rms;
}

void Image::SkyLevel(Pixel *Mean, Pixel *Sigma) const
{
   Frame frame(*this);
   SkyLevel(frame, Mean, Sigma);
}


static double sq(const double &x) { return x*x;}

static void clipped_average_rms(double *pixelValues, const int count,
				double &average, double &rms)
{
  double median = DArrayMedian(pixelValues, count);
#ifdef STORAGE
  double *squares = new double [count];
  
  for (int i=0; i < count ; ++i)
    squares[i] = sq(pixelValues[i]-median);
  
  double var = DArrayMedian(squares,count);
  if (var > 0)  rms = sqrt(var); else rms = 0; 
  
  return;
#endif

  average = median;
  rms = 100000.0;
  
  double meanval,sigval;
  int nloop = 10;
  double quantization = 1.;
  for(int loop=0; loop<nloop; loop++)
    {
      meanval = 0.0;
      sigval = 0.0;
      int nval = 0;
      for (int i=0; i < count; i++)
	{
	  double val = pixelValues[i] - average;
	  if (i!=0) 
	    // compute difference to previous value (they are sorted) to evaluate data quantization.
	    {
	      double diff = pixelValues[i] - pixelValues[i-1];
	      if (diff != 0) quantization = min(diff,quantization);
	    }
	  if (fabs(val)<4*rms) 
	    { // Fix for images with very small sigma (eventually 0)
	      meanval += val;
	      sigval += val*val;
	      nval++;
	    }
	}
      /* meanval is the difference to average, see just above */
      meanval /= nval;
      average += meanval;
      sigval = (sigval/(nval-1)) - (meanval*meanval);
      if (sigval > 0)
	{
	  sigval = sqrt(sigval);
	  rms = sigval;
	}
      else sigval = sqrt(sigval/(nval-1));
      if (meanval == 0) break; // why?
      if (rms < quantization) break; // no need to loop
    }
}


static int greatest_common_divider(int a, int b)
{
  if (a<b) swap(a,b);
  do {
    int q = a/b;
    int r = a-b*q;
    if (r==0) return b;
    a = b;
    b = r;
  } while ( b!= 1);
  return 1;
}

void Image::SkyLevel(const Frame &AFrame, Pixel *Mean, Pixel *Sigma) const
{
  // is AFrame inside the image ?
  Frame frame = AFrame*Frame(*this);
  int npix = int(frame.Nx()*frame.Ny());
  int nvalues=10000;
  int step = npix/nvalues;
  if (npix < nvalues)
    {
      step = 1;
      nvalues = npix;
    }
  
  // lower step until it is prime with nx, unless we sample columns
  int nx = int(frame.Nx());
  int ny = int(frame.Ny());
  while (greatest_common_divider(nx,step) != 1)
    step--;

  //update nvalues accordingly;
  nvalues = npix/step;

  double *pixelValues = new double[nvalues];
  int count =0;
  int imin = int(frame.xMin);
  int jmin = int(frame.yMin);
  int di=0;
  int dj=0;
  while( count < nvalues)
    {
      int i = imin + di;
      int j = jmin + dj;
      if (frame.InFrame(Point(i,j)))
	{
	  pixelValues[count] = (*this)(i,j);
	  ++count;
	}
      int toto = di+step;
      dj += toto/nx;
      di = toto%nx;
      if (dj>=ny) break;
    }

  double average, rms;
  //  if (count != nvalues) cout << " ca ne va pas du tout " << endl;
  clipped_average_rms(pixelValues, count, average, rms);
   
  *Mean = average; *Sigma = rms;
  delete [] pixelValues;
}


#ifdef STORAGE
void Image::convolve(Kernel& K, Image& Result)
{
int ii, ij;
int ik,jk;
int skx = K.HSizeX();
int sky = K.HSizeY();
Pixel sum;

for (ii=sky; ii<=ny-sky; ii++)
  for (ij = skx; ij <nx-skx; ij++)
    {
    sum = 0;
    for (ik = -sky; ik <= sky; ik++)
    for (jk = -skx; jk <= skx; jk++)
      sum += K(jk,ik)*value(ij-jk,ii-ik);
    Result.value(ij,ii) = sum;
    }
}


void Image::fill(FlFunc func)
{int i,j;
for (i=minx(); i <= maxx();i++) for (j=miny(); j <= maxy(); j++)
  value(i,j) = func(i,j);
}
#endif

#include <stdio.h>

void Image::dump() const 
{
int i,j;
for (i=minx(); i<= maxx(); i++)
  {
  for (j=miny(); j<= maxy(); j++)    printf("%5.3f ",(*this)(i,j));
  printf("\n");
  }
printf("\n");
}


static Pixel hmedian(Pixel *array, int size)
{
sort(array, array+size);
return size&1? array[size/2] : (array[size/2-1] + array[size/2])/2.0;
}

void Image::MedianFilter(const int HalfWidth)
{
  int span = 2*HalfWidth+1;
  int npix = span*span;
  Pixel *array = new Pixel[npix];
  // cout << "before median filtering : min" << MinValue() << " max " << MaxValue() << endl;
  Image filtered(nx,ny);
  for (int j=0; j<ny; ++j)
    {
      for (int i=0; i<nx; ++i)
	{
	  npix = 0;
	  for (int jj = max(0,j-HalfWidth); jj <= j+HalfWidth; ++jj)
	    {
	      if (jj == ny) break;
	      for (int ii= max(0, i-HalfWidth); ii <= i+HalfWidth; ++ii)
		{
		  if (ii==nx) break;
		  array[npix++] = value(ii,jj);
		}
	    }
	  filtered(i,j) = hmedian(array,npix);
	}
    }
  delete [] array;
  *this = filtered;
  // cout << "after median filtering : min" << MinValue() << " max " << MaxValue() << endl;
}


void Image::Surface(const int MeshStep,Image & Result)
{
  Pixel mean, sigma ;
  // Le sigma est mal calcule avec MeanSigValue
  MeanSigmaValue(&mean, &sigma);
  cout << "Mean, Sigma : " << mean << " " << sigma << endl ;

  ImageBack back(*this,MeshStep);
  for (int j = 0; j< ny; j++)
    for (int i = 0; i < nx; i++) 
      {
	Result.value(i,j) = value(i,j) - back.Value(i,j);
      }  

  Result.MeanSigmaValue(&mean, &sigma);
  cout << "Mean, Sigma : " << mean << " " << sigma << endl ;
}

double EmpiricCov(const Image &Im1,const Image &Im2)
  // Returns the empirical covariance of 2 images.
{
  if (Im1.Nx() != Im2.Nx() || Im1.Ny() != Im2.Ny() ) 
    {
      cerr << "EmpiricCov : Image 1 and 2 are not the same size. Returning 0" << endl;
      return 0;
    }
  double sum1 = 0;
  double sum2 = 0;
  double sum12 = 0;
  long npix = Im1.Nx()*Im1.Ny();
  Pixel *p1, *p2;
  Pixel *end = Im1.end();
  for (p1=Im1.begin(), p2=Im2.begin(); p1<end; ++p1,++p2)
    {
      sum1 += *p1;
      sum2 += *p2;
      sum12 += (*p1)*(*p2);
    }
  sum1 /= npix;
  sum2 /= npix;
  sum12 /= npix;
  return (sum12 - sum1*sum2);
}



static Pixel local_median(const Image &Img, const int ic, const int jc, 
			  const int size )
{
  int npix = size*size;
  Pixel *buf = new Pixel[npix];
  int count = 0;
  
  for (int j = jc-1; j <= jc+1; ++j)
    for (int i = ic-1; i <= ic+1; ++i)
      buf[count++] = Img(i,j);
  
  Pixel med =  FArrayMedian(buf, npix);
  delete [] buf;
  return med;
}

//Laplacian filter
/*!Cuts (based on the article -> astro-ph/0108003): 
  -cut_lap : the laplacian operator increases the noise by a factor of 
  "sqrt(1.25)"
  
  -cut_f : 2*sigma(med), where sigma(med) is the variance of the
  sky's median calculated in a box (3*3), here. 
  (sigma(med) = sigma(sky)*1.22/sqrt(n); n = size of the box)
  
  -cut_lf : calculated from the article.
  Factor 2.35 -> to have the seeing in arc sec */

int Image::LaplacianFilter(const double & Sigma, const double &Mean,
			   const double &seeing, Image &CosmicImage)
{
  int xmax = this->Nx()-1;
  int ymax = this->Ny()-1;
  float l, f, med;
  
  float cut = 4 * Sigma;
  float cut_lap = cut * sqrt(1.25);
  float cut_f = 2 * (Sigma*1.22/3);
  float cut_lf = 2./(seeing*2.35-1);
  
  int count = 0;
  
  for (int j = 1; j < ymax; j++)
    for (int i = 1; i < xmax; i++)
      {
	//Calculation of the laplacian and the median only for pixels > 3 sigma
	if ((*this)(i,j) > cut) 
	  {
	    l = (*this)(i,j) - 0.25*((*this)(i-1,j) + (*this)(i+1,j) 
				     + (*this)(i,j-1) + (*this)(i,j+1));
	    med = local_median(*this, i, j, 3);
	    f = med - Mean;  //f is invariant by addition of a constant
	    
	    //Construction of a cosmic image
	    if( l>cut_lap && ( f<cut_f || (l/f)>cut_lf) )
	      {
		CosmicImage(i,j) = 1;
		(*this)(i,j) = med;
		count++;
	      }
	    // Image is set to 0 by default
	  }
      }
  CosmicImage.Simplify(0.5);
  return count;
}

void Image::Cosmics(const double &Sigma, const double &Mean,
		    const double &seeing, Image &CosmicImage)  
{
  int iter = 0, count = 1000000;
  
  while (count && iter <5) 
    {
      cout << " Iter " << iter+1;
      count = LaplacianFilter(Sigma, Mean, seeing, CosmicImage);
      cout << " Number of cosmic found " << count << endl;
      iter++;
    }
}

double ImageAndWeightError(const Image &I, const Image &W, const double MinWeight, double *AverageShift)
{
  same_sizes(I,W);
  int npixTot = I.Nx()*I.Ny();
  int npix = min(10000, npixTot);
  int step = max(npixTot/npix, 1);
  double *pixels = new double[npix];
  Pixel *pi = I.begin();
  Pixel *pw = W.begin();
  int count=0;
  for (int i=0; i < npix; ++i)
    {
      if (*pw > MinWeight)
	pixels[count++] = *pi * sqrt(*pw);// only keep significant pixels
      pi += step;
      pw += step;
    }
  npix = count;
  double average,rms;
  clipped_average_rms(pixels, npix, average, rms);
  delete [] pixels;
  if (AverageShift) *AverageShift = average;
  return rms;
}


