#include "image.h"
#include "gtransfo.h"
#include "frame.h"
#include "imageutils.h"
#include "imageinterpolation.h"

Image IntegerShiftImage(const Image& inputimage, const Gtransfo & g, int nx, int ny, float DefaultVal)
{
  double endXStart = double(inputimage.Nx()-1);
  double endYStart = double(inputimage.Ny()-1);
  
  if (nx==0) nx=inputimage.Nx();
  if (ny==0) ny=inputimage.Ny();  
  Image result(nx,ny);

  for (int j=0; j<ny; j++) 
    for (int i=0; i<nx; i++)
      {
	double xout,yout;
	g.apply(i, j, xout, yout);
	if (xout >= 0.0 && xout < endXStart && yout >= 0.0 && yout < endYStart)
	  result(i,j) = inputimage(int(floor(xout+0.5)), int(floor(yout+0.5)));
	else
	  result(i,j) = DefaultVal;

      }
  
  return result;
}

Image GtransfoImage(const Image& inputimage, const Gtransfo & g, int nx, int ny, 
			   float DefaultVal, const int interpLevel,
			   const bool IsVarianceMap)
{
  
  if (IsIdentity(&g)) return inputimage;
  if (IsIntegerShift(&g)){
    cout << " GtransfoImage: transfo is an integer shift, switching to no interpolation\n";
    return IntegerShiftImage(inputimage, g, nx, ny, DefaultVal);
  }

  double endXStart = double(inputimage.Nx()-1);
  double endYStart = double(inputimage.Ny()-1);
  
  if (nx==0) nx=inputimage.Nx();
  if (ny==0) ny=inputimage.Ny();
  
  Image result(nx,ny);
  
  double averageJacobian = fabs(g.Jacobian(int(nx/2),int(ny/2)));

  for (int j=0; j<ny; j++) 
    for (int i=0; i<nx; i++)
      {
	double xout,yout;
	g.apply(i,j,xout,yout);
	if (xout >= 0.0 && xout < endXStart && yout >= 0.0 && yout < endYStart)
	  result(i,j) = Interpolate(inputimage,xout,yout,interpLevel,IsVarianceMap) 
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

// still used as a check, as saturation value are still sometimes misset.
#include "histo1d.h"

double ComputeSaturation(const Image& image)
{
  double maxval = image.MaxValue();
  Histo1d histo(int(maxval/500.), 0, maxval + 1); //HC
  for (int j = 0; j < image.Ny() ; ++j)
    for (int i = 0; i < image.Nx() ; ++i)
      histo.Fill(image(i,j), 1 );
 
  double X;
  double Y = histo.MaxBin(X);
  double maxcoup = X;
  double Saturation=maxval, Sat;
  int count =1, nombre =0;
  for (int l=0; l<20; l++)
    {
      histo.ZeroBins(maxval*0.05*l, maxval*0.05*(l+1));
      Sat = X;
      Y = histo.MaxBin(X);
      if (Sat == X)
        {
          count ++;
          if (count > nombre && Sat != maxcoup)
            {
              nombre = count;
              Saturation =Sat;
            }
        }
      else
        {
          count =1;
 
        }
    }
  return Saturation;
}


/*!  convolve a mask image consisting of "patches" filled
with a uniform value separated by 0's. we enlarge the patches
with the same value, and put -1 in case of conflicting values */
 
bool ConvolveSegMask(const Image &In, Image &Out, const int ExtraSize, int Ny)
{
  if (!In.SameSize(Out))
    {
      cout << "EnlargeMask : In and Out images should have the same size " << endl;
      return false;
    }
  if (ExtraSize == 0)
    {
      Image *out = &Out;
      *out = In;
      return true;
    }

  Image *localOut = &Out;
  // if In and out are the same image, the algorithm just does not work...
  if (&In == &Out)
    {
      localOut = new Image(Out.Nx(), Out.Ny());
    }

  int nx = In.Nx();
  int ny = In.Ny();
  int b = ExtraSize; // too long a name!
  int span = 2*b+1;
  if (Ny > 0 && Ny < ny)
    ny = Ny ;
  for (int j=b; j<ny-b;++j)
    {
      const Pixel *pi = &In(b,j);
      for (int i=b; i<nx-b;++i,++pi)
        {
          if (*pi == 0) continue;
          Pixel val = *pi;
          for (int jj=j-b; jj <= j+b; ++jj)
            {
	      Pixel *po = &(*localOut)(i-b,jj);
                for (int kk = span; kk; kk--, po++)
                  {
                    if (*po == val) continue;
                    if (*po == 0) {*po = val; continue;}
                    *po = -1;
                  }
            }
        }
     }
  if (&In == &Out)
    {
      Out = *localOut;
      delete localOut;
    }
  return true;
}

  


