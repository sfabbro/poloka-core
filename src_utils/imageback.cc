#include <math.h>
#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif

#include "frame.h"
#include "imageback.h"
#include "imageinterpolation.h"

void  Apply_Median_Filter(Image & backValue, Image & backRms, int filter_half_width) ;

ImageBack::ImageBack(Image const &SourceImage, int MeshStep, 
		     const Image*Weight, int filter_half_width) :
  filterHalfWidth(filter_half_width),
     meshStepX(MeshStep) , 
     meshStepY(MeshStep) , 
     //     sourceImage(SourceImage),
     nxImage(SourceImage.Nx()),
     nyImage(SourceImage.Ny()),
     nx(int(ceil(double(SourceImage.Nx())/double(meshStepX)))),
     ny(int(ceil(double(SourceImage.Ny())/double(meshStepY)))),
     backValue(nx,ny),
     backRms(nx,ny) //,
     //     weight(Weight)
{
  do_it(SourceImage, Weight);
}

ImageBack::ImageBack(Image const &SourceImage, int MeshStepX, 
		     int MeshStepY, const Image* Weight, int filter_half_width) : 
  filterHalfWidth(filter_half_width),
     meshStepX(MeshStepX) , 
     meshStepY(MeshStepY) , 
     //     sourceImage(SourceImage),
     nxImage(SourceImage.Nx()),
     nyImage(SourceImage.Ny()),
     nx(int(ceil(double(SourceImage.Nx())/double(meshStepX)))),
     ny(int(ceil(double(SourceImage.Ny())/double(meshStepY)))),
     backValue(nx,ny),
     backRms(nx,ny) //,
     //     weight(Weight)
{
  do_it(SourceImage,Weight);
}


#include "fitsslice.h"
ImageBack::ImageBack(const string FileSourceName, int MeshStepX, int MeshStepY,
	       const string FileWeightName, int filter_hw) :
  filterHalfWidth(filter_hw), 
     meshStepX(MeshStepX) , 
     meshStepY(MeshStepY) ,
     nxImage(FitsHeader(FileSourceName).KeyVal("NAXIS1")),
     nyImage(FitsHeader(FileSourceName).KeyVal("NAXIS2")),
     nx(int(ceil(double(FitsHeader(FileSourceName).KeyVal("NAXIS1"))/double(meshStepX)))),
     ny(int(ceil(double(FitsHeader(FileSourceName).KeyVal("NAXIS2"))/double(meshStepY)))),
     backValue(nx,ny),
     backRms(nx,ny) 
{ 
 
  FitsSlice SourceImage(FileSourceName,MeshStepY,0);
  FitsSlice WeightImage(FileWeightName,MeshStepY,0);

  do_it_slices(SourceImage,WeightImage);
}



/* most of the algorithm is borrowed from sextractor : http://terapix.iap.fr/sextractor */

#include "frame.h"


static void truncated_mean_sig(const Frame &Pad, const Image &im, 
			       const Image *Weight, 
			       double &mean, double &sigma, int &npix)
{
double sum = 0;
double sum2 = 0;
double sumw = 0;
 int imin = int(Pad.xMin);
 int imax = int(Pad.xMax);
 int jmin = int(Pad.yMin);
 int jmax = int(Pad.yMax);
 double value,w;
 if (Weight)
   {
     for (int j=jmin; j<jmax; ++j)
       for (int i=imin; i < imax; ++i)
	 {
	   value = im(i,j);
	   w = (*Weight)(i,j);
	   sum += w*value;
	   sum2 += w*value*value;
	   sumw += w;
	 }
   }
 else
   {
     for (int j=jmin; j<jmax; ++j)
       for (int i=imin; i < imax; ++i)
	 {
	   value = im(i,j);
	   sum += value;
	   sum2 += value*value;
	   sumw += 1;
	 }
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
     w = 1;
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

   mea = med = 0;
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
      mea += (pix*i);  // should be i+0.5.  corrected downwards in the definition of qzero
      sig +=  (pix*i*i);
      }

    //handling of empty bins (happens on quantized data)
    while (*hilow == 0 && hilow < histo+hcut) hilow++;
    while (*hihigh == 0 && hihigh > histo+lcut) hihigh--;

    med = (hihigh-histo)+
      (hilow-hihigh)*(0.5+(highsum-lowsum)/(2.0*(*hilow>*hihigh?
						 *hilow:*hihigh)));

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


  // mea and med are expressed in channel number where the first channel is labeled 0.
  // to convert to actual abcissa, one has to account for a half bin offset:

   double binSize = Histo.BinWidth();
   double qzero = Histo.Minx()+binSize*0.5;

   // does (fabs(sigma/(sig*binSize)-1) < 0.0) ever happen ??
   //  mean = fabs(sig)>0.0? (fabs(sigma/(sig*binSize)-1) < 0.0 ?
   //            qzero+mea*binSize
   //            :(fabs((mea-med)/sig)< 0.3 ?
   //              qzero+(2.5*med-1.5*mea)*binSize
   //            :qzero+med*binSize))
   //                  :qzero+mea*binSize;                   :qzero+mea*binSize;

   // DH : M. BETOULE PATCH
   if (fabs(sig)>0.0){
     if (fabs((mea-med)/sig)< 0.3) // non crowded
       {
	 mean = qzero+mea*binSize;
       }
     else //crowded
       {
	 mean = qzero+(2.5*med-1.5*mea)*binSize; // Sextractor's approximation of the mode 
       }
   }
   else
     {
       mean = qzero+mea*binSize;
     }

   sigma = sig*binSize;
  }





/* these defines I do not understand (all) are from sextractor */
#define N_SIGMAS 5.  
#define Q_AMIN 4.
#define MAX_BINS 4096
#define MIN_PIX_FRAC 0.1


static void do_one_pad(const Image &SourceImage, const Image *Weight, const Frame &Pad, double &mean, 
				  double &sigma)
{
  int npix;
  truncated_mean_sig(Pad, SourceImage, Weight, mean, sigma, npix);
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
  // DEBUG
  // cout << " nbins,minh,maxh " << histo.Nx() << ' ' << histo.Minx() << ' ' <<  histo.Minx()+histo.Nx()*histo.BinWidth() << endl;

  int imin = int(Pad.xMin);
  int imax = int(Pad.xMax);
  int jmin = int(Pad.yMin);
  int jmax = int(Pad.yMax);
  if (!Weight)
   {
     for (int j = jmin; j< jmax; ++j) for (int i=imin; i< imax; ++i)
       histo.Fill(SourceImage(i,j));
   }
  else
    {
      for (int j = jmin; j< jmax; ++j) for (int i=imin; i< imax; ++i)
	histo.Fill(SourceImage(i,j),(*Weight)(i,j));
    }
  backguess(histo,mean,sigma);
}


void ImageBack::do_it(const Image &SourceImage, const Image* Weight)
{
for (int j=0; j<ny; j++)
  {
  int height = (j<ny-1) ? meshStepY : SourceImage.Ny() - j *meshStepY;
  for (int i=0; i<nx; i++)
    {
    int width = (i<nx-1)? meshStepX : SourceImage.Nx()- i*meshStepX;
    int startx = i*meshStepX;
    double mean, sigma;
    do_one_pad(SourceImage, Weight, Frame(startx, j*meshStepY, 
		     startx+width, j*meshStepY+height),mean, sigma);
    backValue(i,j) = mean;
    backRms(i,j) = sigma;
    }
  }

 Apply_Median_Filter(backValue, backRms, filterHalfWidth) ;
} 








void ImageBack::do_it_slices(FitsSlice & SourceImage, FitsSlice & WeightImage)
{
  // normalement,autant de slices que de ny
  int j = 0 ;
  do {
  
    int height = (j<ny-1) ? meshStepY : (nyImage - j *meshStepY);
    if ( height != SourceImage.SliceSize() )
      {
	cerr << "erreur in slicing at line : " << SourceImage.ImageJ(0) <<  " computed height : " << height << " (j="<< j << " /ny= " << ny << ") slice height " << SourceImage.SliceSize() << endl ;
      }
    for (int i=0; i<nx; i++)
      {
	int width = (i<nx-1)? meshStepX : (SourceImage.Nx()- i*meshStepX);
	int startx = i*meshStepX;
	double mean, sigma;
	do_one_pad(SourceImage, &WeightImage, Frame(startx, 0, startx+width, height),mean, sigma);
	backValue(i,j) = mean;
	backRms(i,j) = sigma;
      }
    j++;
  }
  while (SourceImage.LoadNextSlice() && WeightImage.LoadNextSlice()) ;
  cerr << "Nslices prevues : " << j << " ny : " << ny << endl ;

//
 Apply_Median_Filter(backValue, backRms, filterHalfWidth) ;
}


void  Apply_Median_Filter(Image & backValue, Image & backRms, int filter_half_width)
{
  int nx = backValue.Nx();
  int ny = backValue.Ny();
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
  backValue.MedianFilter(filter_half_width);
  backRms.MedianFilter(filter_half_width);
}


































Image*  ImageBack::BackgroundImage(const Frame &AFrame) const
{
  int iMin = int(floor(AFrame.xMin+0.5));
  int jMin = int(floor(AFrame.yMin+0.5));
  // xmax amd ymax are within the image (see Frame::Frame(const Image&)) 
  int iMax = int(floor(AFrame.xMax+0.5));
  int jMax = int(floor(AFrame.yMax+0.5));
  Image *result = new Image(iMax-iMin+1,jMax-jMin+1);
  for(int j = jMin ; j<=jMax ; j++)
    for(int i = iMin ; i<=iMax ; i++)
      (*result)(i-iMin,j-jMin) = Value(i,j);
  return(result) ;
}

Image*  ImageBack::BackgroundImage() const
{
  return BackgroundImage(Frame(0,0,nxImage-1, nyImage-1));
}




Pixel ImageBack::Value(const int x, const int y) const
{
double xBack = ((x + 0.5)/meshStepX) - 0.5;
double yBack = ((y + 0.5)/meshStepY) - 0.5;

return Interpolate(backValue,xBack,yBack,2);
}

Pixel ImageBack::Rms(const int x, const int y) const
{
double xBack = ((x + 0.5)/meshStepX) - 0.5;
double yBack = ((y + 0.5)/meshStepY) - 0.5;

return Interpolate(backRms,xBack,yBack,2);
}


//By July 2006, these IO functionalities are not used!!
/* 
************************ WARNING *****************************
 if you are considering using them, you should
- have 2 separate routines to write the BackValue and BackRms
- put the nxImage and nyImage in the fits headers
- avoid mangling the user provided filenames
*/


#include "fitsimage.h"
#include <stdio.h>
#include "fileutils.h"

ImageBack::ImageBack(char *FileName) : backValue(), backRms()
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

void ImageSurface(const Image& inputimage, const int MeshStep,Image & Result) {
  Pixel mean, sigma ;
  // Le sigma est mal calcule avec MeanSigValue
  inputimage.MeanSigmaValue(&mean, &sigma);
  cout << "Mean, Sigma : " << mean << " " << sigma << endl ;

  ImageBack back(inputimage,MeshStep);
  for (int j = 0; j< inputimage.Ny(); j++)
    for (int i = 0; i < inputimage.Nx(); i++) 
      {
	Result(i,j) = inputimage(i,j) - back.Value(i,j);
      }  

  Result.MeanSigmaValue(&mean, &sigma);
  cout << "Mean, Sigma : " << mean << " " << sigma << endl ;




}
