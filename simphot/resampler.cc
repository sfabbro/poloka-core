#include "resampler.h"

#include "gtransfo.h"
#include "array4d.h"
#include "cmath" // for floor


/* Before modifying routines here, read the initial comment in resampler.h */


int ResamplerBoundarySize()
{
  return 1;
}


/*! Computes the resampling coefficients. Transfo goes from target to original
(current to model in the framework of SimPhotFit), and Coeffs contains 
the boundaries in the (normally "reduced", i.e. centered on the object) 
coordinates of the current image */
void ResamplerComputeDerivatives(const Gtransfo* Transfo, 
				 Array4D & Coeffs)
{
  double xt,yt, xi, yi;
  int ic,jc;
  // non oversampling case :
  double dx,dy;
  double xf0, xf1, xf2, yf0, yf1, yf2;

  //TODO : account for the jacobian of the Transfo.

  for (int b = Coeffs.ymin; b < Coeffs.ymax; ++b)
    for (int a = Coeffs.xmin; a < Coeffs.xmax; ++a)
      {
	Transfo->apply(double (a), double(b), xt, yt);
	xi = floor(xt+0.5);
	yi = floor(yt+0.5);
	ic = int(xi);
	jc = int(yi);
	CoeffBlock &block = Coeffs(a,b);
	block.Allocate(ic-1,jc-1,ic+2, jc+2);
	dx = xt -xi;
	dy = yt -yi;
   
	xf0 = 0.25*(dx-1)*dx;
	xf1 = -0.5*(dx-1)*(dx+1);
	xf2 = 0.25*(dx+1)*dx;
	yf0 = (dy-1)*dy;
	yf1 = -2*(dy-1)*(dy+1);
	yf2 = (dy+1)*dy;
	// totaly unrolled loops.. should be faster than actual loops.
	// next optimization : push the "write" pointer by hand.
	block(ic-1,jc-1) = xf0*yf0;
	block(ic,  jc-1) = xf1*yf0;
	block(ic+1,jc-1) = xf2*yf0;

	block(ic-1,jc) =   xf0*yf1;
	block(ic,jc) =     xf1*yf1;
	block(ic+1,jc) =   xf2*yf1;

	block(ic-1,jc+1) = xf0*yf2;
	block(ic,jc+1) =   xf1*yf2;
	block(ic+1,jc+1) = xf2*yf2;
      }
}

#include "pixelblock.h"
void ResampleImage(const PixelBlock &In, const Gtransfo* Out2In,
		   PixelBlock &Out)
{
  double x,y,xi,yi,xt,yt;
  int ic,jc;
  double xf0,xf1,xf2,yf0,yf1,yf2;
  double dx,dy;

  for (int j = Out.ymin; j < Out.ymax; ++j)
  for (int i = Out.xmin; i < Out.xmax; ++i)
    {
      x = i;
      y = j;
      Out2In->apply(x,y,xt,yt);
      xi = floor(xt+0.5);
      yi = floor(yt+0.5);
      ic = int(xi);
      jc = int(yi);
      if (ic-In.xmin<1 || In.xmax-ic<2 || jc - In.ymin <1 || In.ymax - jc< 2)
	{
	  Out(i,j) = 0;
	  continue;
	}
      dx = xt-xi;
      dy = yt-yi;
      xf0 = (dx-1)*dx;
      xf1 = (dx-1)*(dx+1);
      xf2 = (dx+1)*dx;
      yf0 = (dy-1)*dy;
      yf1 = (dy-1)*(dy+1);
      yf2 = (dy+1)*dy;
      Out(i,j) = 0.25*(
         (    In(ic-1, jc-1)  *xf0
         -2.0*In(ic  , jc-1)*xf1
           +  In(ic+1, jc-1)*xf2) * yf0
        -2.0*(In(ic-1, jc  )*xf0
         -2.0*In(ic  , jc  )*xf1
            + In(ic+1, jc  )*xf2) * yf1
         +   (In(ic-1, jc+1)*xf0
         -2.0*In(ic  , jc+1)*xf1
            + In(ic+1, jc+1)*xf2) * yf2
	 );
    }
}


void ConvolveImage(const PixelBlock &In, const PixelBlock &Kernel,
		   PixelBlock &Out)
{
  Out.Allocate(In.xmin+Kernel.xmax-1,In.ymin+Kernel.ymax-1,
	       In.xmax+Kernel.xmin, In.ymax+Kernel.ymin);
  for (int d = Out.ymin; d < Out.ymax; ++d)
    for (int c = Out.xmin; c < Out.xmax; ++c)
      {
	double val = 0;
	for (int j= Kernel.ymin; j <Kernel.ymax; ++j)
	  for (int i = Kernel.xmin ; i <Kernel.xmax ; ++i)
	    {
	      val += In(c-i,d-j)*Kernel(i,j);
	      //	      cout << c << ' ' << d << ' ' << i << ' ' << j << ' ' << val << endl;
	    }
	Out(c,d) = val;
      }
}
