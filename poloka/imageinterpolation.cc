#include <cmath>

#include <poloka/image.h>
#include <poloka/imageinterpolation.h>

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

Pixel Interpolate(const Image& inputimage, const double x, const double y, const int level,
			 const bool IsVarianceMap)
{
  int xp,yp;
  int nx = inputimage.Nx();
  int ny = inputimage.Ny();
  
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

      if (nx == 1)
	{
	  if (ny == 1) return inputimage(xp,yp);
	  else
	    return inputimage(xp,yp)*(1-dy) + inputimage(xp,yp+1)*dy;
	}
      if (ny == 1) // then nx != 1 
	return (inputimage(xp,yp)*(1-dx) + inputimage(xp+1,yp)*dx);

      if (!IsVarianceMap)
	{
	  return ((inputimage(xp,yp)*(1-dx) + inputimage(xp+1,yp)*dx)*(1-dy) +
		  (inputimage(xp,yp+1)*(1-dx) + inputimage(xp+1,yp+1)*dx)*dy);
	}
      else // transforming variances,
	{ // the same as above with squared coefficients:
	  return (
	   (inputimage(xp,yp)*(1-dx)*(1-dx) + 
	    inputimage(xp+1,yp)*dx*dx)
	   *(1-dy)*(1-dy) 
	   +(inputimage(xp,yp+1)*(1-dx)*(1-dx) + 
	     inputimage(xp+1,yp+1)*dx*dx)
	   *dy*dy);
	}
    }
  else // Quadratic interpolate
    {
      return QuadraticInterpolate(inputimage, x, y, xp, yp, IsVarianceMap);
    }
}

