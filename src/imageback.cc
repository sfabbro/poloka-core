#include <math.h>
#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif

#include "imageback.h"

ImageBack::ImageBack(Image const &SourceImage, int MeshStep, 
		     const Image*Weight) : 
     meshStepX(MeshStep) , 
     meshStepY(MeshStep) , 
     sourceImage(SourceImage),
     nx(int(ceil(double(sourceImage.Nx())/double(meshStepX)))),
     ny(int(ceil(double(sourceImage.Ny())/double(meshStepY)))),
     backValue(nx,ny),
     backRms(nx,ny),
     weight(Weight)
{
do_it();
}

ImageBack::ImageBack(Image const &SourceImage, int MeshStepX, 
		     int MeshStepY, const Image* Weight) : 
     meshStepX(MeshStepX) , 
     meshStepY(MeshStepY) , 
     sourceImage(SourceImage),
     nx(int(ceil(double(sourceImage.Nx())/double(meshStepX)))),
     ny(int(ceil(double(sourceImage.Ny())/double(meshStepY)))),
     backValue(nx,ny),
     backRms(nx,ny),
     weight(Weight)
{
do_it();
}





/* most of the algorithm is borrowed for sextractor : http://terapix.iap.fr/sextractor */

#include "frame.h"


static void truncated_mean_sig(const Frame &Pad, const Image &im, 
			       const Image *Weight, 
			       double &mean, double &sigma, int &npix)
{
double sum = 0;
double sum2 = 0;
double sumw = 0;
double value,w;
 int imin = int(Pad.xMin);
 int imax = int(Pad.xMax);
 int jmin = int(Pad.yMin);
 int jmax = int(Pad.yMax);
 for (int j=jmin; j<jmax; ++j)
   for (int i=imin; i < imax; ++i)
     {
       value = im(i,j);
       w = Weight? (*Weight)(i,j) : 1;
       sum += w*value;
       sum2 += w*value*value;
       sumw += w;
  }
 double sig;
 if (sumw)
   {
     mean  = sum/sumw;
     sig = (sum2/sumw - mean*mean); if (sig>0) sig = sqrt(sig); else sig = 0;
   }
 else { sigma = 1; mean = 0; npix = 10; return;}
 int loop;
 for (loop=2; loop>0; loop--)
   {
     double lcut =  -2.0*sig;
     double hcut =  2.0*sig;
     sum = sum2 = 0; sumw = 0; npix = 0;
     for (int j=jmin; j<jmax; ++j)
       for (int i=imin; i < imax; ++i)
	 {
	   value = im(i,j)-mean;
	   if (Weight) w = (*Weight)(i,j);
	   if (w == 0) continue;
	   if (value > hcut || value < lcut) continue;
	   sum += value*w;
	   sum2 += value*value*w;
	   npix++;
	   sumw += w;
	 }
     if (sumw)
       {
	 sum =  sum/sumw;
	 sum2 = sum2/sumw - sum*sum;
	 mean += sum;
	 sig = (sum2>0) ? sqrt(sum2) : 0.0;
       }
     else break; // keep "old" values
   }
 sigma = sig;
}


#define BIG 1e30
#define BAD -(BIG)

#include "histo1d.h"

static void  backguess(const Histo1d &Histo, double &mean, double &sigma)
#define EPS (1e-4)  /* a small number */

  {
   const float    *histo = Histo.array();
   const float     *hilow, *hihigh, *histot;
   double   lowsum, highsum, sum;
   double   ftemp, mea, sig, sig1, med;
   int      i, n, lcut,hcut, nlevelsm1;


  hcut = nlevelsm1 = Histo.Nx() - 1;
  lcut = 0;

  sig = 10.0*nlevelsm1;
  sig1 = 1.0;
  for (n=100; n-- && (sig>=0.1) && (fabs(sig/sig1-1.0)>EPS);)
    {
    sig1 = sig;
    sum = mea = sig = 0.0;
    lowsum = highsum = 0;
    histot = hilow = histo+lcut;
    hihigh = histo+hcut;
    for (i=lcut; i<=hcut; i++)
      {
      if (lowsum<highsum)
        lowsum += *(hilow++);
      else 
        highsum +=  *(hihigh--);
      double pix = *(histot++);
      sum += pix;
      mea += (pix*i);
      sig +=  (pix*i*i);
      }

    med = (hihigh-histo)+0.5+((double)highsum-lowsum)/(2.0*(*hilow>*hihigh?
                        *hilow:*hihigh));
    if (sum)
      {
      mea /= (double)sum;
      sig = sig/sum - mea*mea;
      }
    sig = sig>0.0?sqrt(sig):0.0;
    lcut = (ftemp=med-3.0*sig)>0.0 ?(int)(ftemp>0.0?ftemp+0.5:ftemp-0.5):0;
    hcut = (ftemp=med+3.0*sig)<nlevelsm1 ?(int)(ftemp>0.0?ftemp+0.5:ftemp-0.5)
                                : nlevelsm1;
    }

   double binSize = Histo.BinWidth();
   double qzero = Histo.Minx();
   mean = fabs(sig)>0.0? (fabs(sigma/(sig*binSize)-1) < 0.0 ?
               qzero+mea*binSize
                :(fabs((mea-med)/sig)< 0.3 ?
                  qzero+(2.5*med-1.5*mea)*binSize
                 :qzero+med*binSize))
                       :qzero+mea*binSize;

   sigma = sig*binSize;
  }





/* these defines I do not understand (all) are from sextractor */
#define N_SIGMAS 5.  
#define Q_AMIN 4.
#define MAX_BINS 4096
#define MIN_PIX_FRAC 0.1


void ImageBack::do_one_pad(const Frame &Pad, double &mean, 
				  double &sigma)
{
  int npix;
  truncated_mean_sig(Pad, sourceImage, weight, mean, sigma, npix);
  if (npix < MIN_PIX_FRAC * Pad.Area())
    {
      mean = BAD;
      sigma = BAD;
      return;
    }
  double step = sqrt(2./M_PI)*N_SIGMAS*Q_AMIN;
  int nbins = min(int(step*npix+1.), MAX_BINS);
  if (sigma<=0) sigma = 1.;
  Histo1d histo(nbins,mean - N_SIGMAS*sigma, mean+N_SIGMAS*sigma);
  int imin = int(Pad.xMin);
  int imax = int(Pad.xMax);
  int jmin = int(Pad.yMin);
  int jmax = int(Pad.yMax);
  if (!weight)
   {
     for (int j = jmin; j< jmax; ++j) for (int i=imin; i< imax; ++i)
       histo.Fill(sourceImage(i,j));
   }
  else
    {
      for (int j = jmin; j< jmax; ++j) for (int i=imin; i< imax; ++i)
	histo.Fill(sourceImage(i,j),(*weight)(i,j));
    }
  backguess(histo,mean,sigma);
}


void ImageBack::do_it()
{
for (int j=0; j<ny; j++)
  {
  int height = (j<ny-1) ? meshStepY : sourceImage.Ny() - j *meshStepY;
  for (int i=0; i<nx; i++)
    {
    int width = (i<nx-1)? meshStepX : sourceImage.Nx()- i*meshStepX;
    int startx = i*meshStepX;
    double mean, sigma;
    do_one_pad(Frame(startx, j*meshStepY, 
		     startx+width, j*meshStepY+height),mean, sigma);
    backValue(i,j) = mean;
    backRms(i,j) = sigma;
    }
  }

// fill in bad values (where the weight is 0 on (almost) a whole pad)
 Image back2(backValue); 
  for (int j=0; j<ny; ++j)
    for (int i=0; i<nx; ++i)
      if (backValue(i,j) <= BAD)
        {
/*------ Seek the closest valid mesh */
	  double d2min = BIG;
	  int nmin = 0;
	  double val = 0;
	  double sval = 0;
	  /* could probably do it more efficiently... for the moment
	     traverse the whole back image: it is not that large usually */
	  for (int jgood=0; jgood<ny; ++jgood)
	  for (int igood =0; igood<nx; ++igood)
	    {
	      if (backValue(igood,jgood)> BAD)
		{
		  double d2 = (i-igood)*(i-igood)+(j-jgood)*(j-jgood);
		  if (d2<d2min)
		    {
		      val = backValue(igood,jgood);
		      sval = backRms(igood,jgood);
		      nmin = 1;
		      d2min = d2;
		    }
		  else if (d2==d2min)
		    {
		      val += backValue(igood,jgood);
		      sval += backRms(igood,jgood);
		      nmin++;
		    }
		}
	    }
	  if (nmin !=0)
	    {
	      back2(i,j) = val/nmin;
	      backRms(i,j) = sval/nmin;
	    }
	  else
	    {
	      back2(i,j) = 0;
	      backRms(i,j) = 1.0;
	      cerr << " did not find a correct pad to copy when filling in bad areas of back image " << endl;
	    }

	}
  backValue = back2;

  // sextractor applies a 'median' filter on average and sigma maps (with a half width of one)
  backValue.MedianFilter(1);
  backRms.MedianFilter(1);
} 

Image*  ImageBack::BackgroundImage()
{
  int tx = sourceImage.Nx();
  int ty = sourceImage.Ny();
  Image *result = new Image(tx,ty);
  for(int i = 0 ; i < tx ; i++)
    for(int j = 0 ; j < ty ; j++)
      (*result)(i,j) = Value(i,j);
  return(result) ;
}

Pixel ImageBack::Value(const int x, const int y) const
{
double xBack = ((x + 0.5)/meshStepX) - 0.5;
double yBack = ((y + 0.5)/meshStepY) - 0.5;

return backValue.Interpolate(xBack,yBack,2);
}

Pixel ImageBack::Rms(const int x, const int y) const
{
double xBack = ((x + 0.5)/meshStepX) - 0.5;
double yBack = ((y + 0.5)/meshStepY) - 0.5;

return backRms.Interpolate(xBack,yBack,2);
}


#include "fitsimage.h"
#include <stdio.h>
#include "fileutils.h"

ImageBack::ImageBack(char *FileName, const Image &Source) : sourceImage(Source), backValue(), backRms()
{
char file_name[80];
sprintf(file_name,"%s_mb",FileName);
// if (!FileExists(FileName)) {is_valid = 0; return;}
FitsImage back(file_name);
nx = back.Nx();
ny = back.Ny();
meshStepX = back.KeyVal("MSTEPX");
meshStepY = back.KeyVal("MSTEPY");
backValue = back;
sprintf(file_name,"%s_ms",FileName);
backRms = FitsImage  (file_name);
}

int ImageBack::Write(string in_name)
{
string out_name = in_name+"_mb";
FitsImage back(out_name, backValue);
back.AddKey("MSTEPX",meshStepX);
back.AddKey("MSTEPY",meshStepY);
out_name = in_name + "_ms";
FitsImage rms(out_name, backRms);
return 1;
}
