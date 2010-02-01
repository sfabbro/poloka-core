#include <iostream>
#include <string>
#include <math.h>

#include "toadscards.h"
#include "fitsimage.h"
#include "convolution.h"
#include "detection.h"
#include "frame.h"

static double sqr(double x){return x*x;}

#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif

//#define DEBUG

struct DatDet{
  double nRadAper;
  double nRadWeight;
  double nSigDet;
  double nSig;
  double nRadBad;
  double nMinDistDouble;
  double compFactor;
  double nbackSideMin;
  double nbackSideMax;
  string fileName;

  DatDet(const string &DatacardsFileName); 
  void Print(ostream & s=cout) const ;
};

#include "fileutils.h" // for FileExists
#include "datacards.h"

DatDet::DatDet(const string &DatacardsFileName)
{
  //default val
  nRadAper = 2.5; // i.e. intg. radius of aperture = 1FWHM
  nRadWeight = 1.6 * nRadAper;
  nRadBad = 2.;
  nSigDet = 2.0; // cut in S/N on the convolved image.
  nSig = 2.5; // (final) cut in S/N on finale flux estimate.
  nMinDistDouble = nRadBad;
  compFactor = 3.; 
  nbackSideMin =  nRadWeight*1.1 ;  
  nbackSideMax = 8 ;
  if (FileExists(DatacardsFileName))
    {
      fileName = DatacardsFileName;
      // overwrite if found in datacards
      DataCards cards(DatacardsFileName);
      if (cards.HasKey("DATDET_NRADBAD"))
	nRadBad = cards.DParam("DATDET_NRADBAD");
      if (cards.HasKey("DATDET_NRADWEIGHT"))
	nRadWeight = cards.DParam("DATDET_NRADWEIGHT");
      if (cards.HasKey("DATDET_NSIG_DET"))
	nSigDet = cards.DParam("DATDET_NSIG_DET");
      if (cards.HasKey("DATDET_NSIG"))
	nSig = cards.DParam("DATDET_NSIG");
      if (cards.HasKey("DATDET_NRADAPER"))
	nRadAper = cards.DParam("DATDET_NRADAPER");
      if (cards.HasKey("DATDET_NMINDISTDOUBLE"))
	nMinDistDouble = cards.DParam("DATDET_NMINDISTDOUBLE");
      if (cards.HasKey("DATDET_BACKRINGMIN"))
	nbackSideMin = cards.IParam("DATDET_BACKRINGMIN");
      if (cards.HasKey("DATDET_BACKRINGMAX"))
	nbackSideMax = cards.IParam("DATDET_BACKRINGMAX");
      if (cards.HasKey("DATDET_COMPFACTOR"))
	compFactor = cards.DParam("DATDET_COMPFACTOR");
    }
  else
    {
      cerr << " Cannot read " << DatacardsFileName << endl;
    }
  if ( nbackSideMax == 0)
    nbackSideMin = 0 ;
  else
    {
      if ( nbackSideMax <nbackSideMin ) swap(nbackSideMin,nbackSideMax);
      if ( nbackSideMin < nRadWeight) // min rad max OR min max rad
	{
	  cerr << "Inner Size of annulus to compute back > radius for weighted flux " << endl ;
	  nbackSideMin = nRadWeight; // min=rad max OR max rad=min
	  if ( nbackSideMax < nRadWeight) // max rad=min
	    {
	      cerr << "Outer Size of annulus to compute back < radius for weighted flux, modifying annulus radii values" << endl ;
	      nbackSideMax = 2.*nRadWeight; // rad=min max=2.*rad
	    }
	}
    }
}


void DatDet::Print(ostream & s) const 
{
  s  << "***** DatDet ******" << endl ;
  s  << " datacards file name " << fileName << endl; 
  s  << " Integration radius for aperture flux (in sig seeing) : " 
     <<  nRadAper << endl;
  s  << " Integration radius for weighted flux (in sig seeing): " 
     <<  nRadWeight << endl;
  s  << " Detection cut in S/N "   << nSigDet << endl;
  s  << " (final) Flux cut in S/N "   << nSig << endl;
  s  << " distance cut (object bad pix) in sig seeing " << nRadBad << endl;
  s  << " min distance (in sig seeing) between 2 detections " 
     << nMinDistDouble << endl;

  s  << " inner / outer radii for the back estimation \"annulus\" (in sig seeing)" << nbackSideMin << ", " << nbackSideMax
     << endl;
  s  << " ********** end of DatDet ************ " << endl;
}



// ----------------------------------- 
// ------ utility routines ------
// ----------------------------------- 

static bool IsLocalMax( Image const & img , int x , int y , int half_square_size)
{
  int tx = img.Nx() ;
  int ty = img.Ny() ;
  float val=0;
  if ( ( x > 0 ) && ( x < tx ) &&
       ( y > 0 ) && ( y < ty ) )
    val = img(x,y);
  int i , j ;
  
  for (i=-half_square_size; i<=half_square_size; i++)
    for (j=-half_square_size; j<=half_square_size; j++)
      {
	if ((i==0) && (j==0)) continue;
	int xx = x + i ;
	int yy = y + j ;
	if ( ( xx > 0 ) && ( xx < tx ) &&
	     ( yy > 0 ) && ( yy < ty ) )
	  if (img(xx,yy)>val) 
	    return false;
      }
  return true;
}


/****************** DetectionProcess *********************/

#define PrintDatacards 1
#define PrintCuts 2



DetectionProcess::DetectionProcess(const string &ImageName, 
				   const string &WeightName, 
				   const double Seeing, const double SigFilter)
{
  
#ifdef DEBUG
  cout << "DetectionProcess::DetectionProcess" << endl;
#endif
  seeing = Seeing;
  sigFilter = SigFilter;
  weightName = WeightName;
  imageName = ImageName;

  weight = NULL;
  imagecv = NULL;
  image = new FitsImage(imageName);

  cvImageNormalizedSig = 0;
  cvImageNormalizedMean = 0;
  /* in principle this one is useless since there is a call in all
     public service routines (with proper print instructions) */
  SetCuts(); 
}


void DetectionProcess::SetCuts(const int ToPrint)
{
  
#ifdef DEBUG
  cout << "DetectionProcess::SetCuts" << endl;
#endif

  //collect cuts in datacards
  DatDet datdet(DefaultDatacards());
  if (ToPrint & PrintDatacards) datdet.Print();
  radAper = seeing * datdet.nRadAper;
  radWeight = seeing *datdet.nRadWeight;
  nSig = datdet.nSig;
  nSigDet = datdet.nSigDet;
  radBad = datdet.nRadBad * seeing;
  minDistDouble = datdet.nMinDistDouble * seeing;
  radBackMin = datdet.nbackSideMin * seeing;
  radBackMax = datdet.nbackSideMax * seeing;
  if (ToPrint & PrintCuts)
    {
      cout << " cuts for image " << imageName << endl;
      cout << " hopefully cutting detections at  " 
	   << nSigDet << " then " << nSig << " sigmas " << endl;
      cout << " radAper = " << radAper << endl;
      cout << " radWeight = " << radWeight << endl;
      cout << " radBad = " << radBad << endl;
      cout << " min dist between 2 det " << minDistDouble << endl;
      cout << " annulus inner/outer radii for back estimation = " 
	   << radBackMin << " / " <<  radBackMax << endl;
      cout << " ******************************" << endl;
    }
}

DetectionProcess::~DetectionProcess()
{
  if (image) delete image;
  if (imagecv) delete imagecv;
  if (weight) delete weight;
}


static void print_stat(const Image &Im, const string &Mess)
{
  Pixel mean, sigma;
  Im.SkyLevel(&mean, &sigma);
  Pixel min,max;
  Im.MinMaxValue(&min,&max);
  cout << Mess << "  m,s: "  << mean << ' ' << sigma 
       << " min " << min << " max " << max << endl;
}




//TO DO : test this routine on synthetic noise data, compute analytically
// what the output should be (Irwin 85 may help) and compare to what comes out.
/* principle of this routine:

the optimal flux estimator for a point like stuff on a weighted
image (i.e. an image that has pixel weights associated) reads:

        sum( psf * image * weight) / sum(psf)
flux = --------------------------------------
          sum(psf^2 *weight) / [sum(psf)]^2
             
Var(flux) =  [sum(psf)]^2 / sum(psf^2 *weight)

This is what is computed in SetDetectionScores, and ImagesConvolve
should just compute the same thing for ecvery pixel, assuming that a 
point source lies in this pixel. In fact what we need to actually
detect is flux/sqrt(var) which reads


flux           sum( psf * image * weight) 
------- = -----------------------------  (1)
sig_flux      sqrt( sum(psf^2 * weight) )

This expression is the significance of the flux (would be) excess in a given
pixel. The sums extend to surrounding pixels, and the psf is to be understood
as the psf valu at the summation pixel, centered on the reference pixel (i.e.)
where we are estimating the flux. SO if we distinguish convolutions (%*%)
from products (*), we get:

flux           [psf %*% (image * weight) ](i)
------- (i) = -----------------------------      (2)
sig_flux      sqrt( (psf^2 %*% weight) (i) )




This expression coincides with what one can readily compute for 
equal weight pixels.
This expression is very however different from what I coded in a first attempt:

      sum(psf * image)/ sum(psf^2)
      -----------------
   sqrt(sum(psf^2 / weights)) / [sum(psf^2)]^2

*/

void DetectionProcess::ImagesConvolve()
{
  /* to compute expresion (3) above we need 2 images : 1 for the numerator
     and 1 for the denominator. we could use the space reserved for 
     the image for (e.g.) the numerator and restore the values at the end 
     of this routine. We do not follow this memory saving approach for now.
  */
  imagecv = new Image(*image);
  // read the weights
  Image weightcv;
  { /* take a copy of the image to delete the FitsImage (and close the file)
       immediately after reading */

    FitsImage w(weightName);
    weightcv = w;
  }
  *imagecv *= weightcv;

  Image kernel = ConvoleGauss1D1D(*imagecv, sigFilter, sigFilter, 
				  /*prec = */ 1.e-3);

  // in principe the kernel is normalized, but who nows how software evolves:
  double fact = 1./kernel.SumPixels();
  *imagecv *= fact;

  // done for the numerator.

  //denominator
  Image kernweight = ConvoleGauss1D1D(weightcv, sigFilter/sqrt(2.), 
				    sigFilter/sqrt(2.), 
				    1.e-3);

  /* there is an issue here: the convolution routines normalize the kernel.
     We want :
     psf^2 %*% weights / [sum (psf)]^2
     and we get 
     psf^2 %*% weights / sum (psf ^ 2)
     the correction factor is then
     sum(psf^2) / [ sum(psf) ]^2
     with the caveat that kernweight == psf^2 
  */


  print_stat(weightcv , " weight after conv");

  fact = kernweight.SumPixels();
  kernweight.ApplyFun(sqrt);
  fact /= sqr(kernweight.SumPixels());
  weightcv *= fact;

  print_stat(weightcv , " weight after conv+norm");

  // we now have the numerator and denominator of expression (2)

// compute the frame in which the convolution was actually computed.
  Frame imFrame(*imagecv);
  int xMarginSize = max(kernel.Nx(), kernweight.Nx())/2;
  int yMarginSize = max(kernel.Ny(), kernweight.Ny())/2;
  imFrame.CutMargin(xMarginSize, yMarginSize);

  // compute the imagecv average:
  Pixel meanIm, sigIm;
  imagecv->SkyLevel(imFrame, &meanIm, &sigIm);

  // compute the weightcv average:
  Pixel meanW, sigW;
  weightcv.SkyLevel(imFrame, &meanW, &sigW);

  // default value for imagecv/weightcv when weightcv =0;
  Pixel defaultVal = meanIm/meanW;

  
  // compute expression (2)
  Pixel *pw = weightcv.begin();
  Pixel *pend = imagecv->end();
  for (Pixel *p = imagecv->begin(); p < pend; ++p, ++pw) 
    {
      if (*p == 0) continue;
      if (*pw == 0) { *p = defaultVal; continue;}
      (*p) /= (*pw);
    }

  // set pixels wher the convolution did not go to defaultVal
  imagecv->Masking(imFrame, defaultVal);

  //compute noise
   cvImageNormalizedSig = ImageAndWeightError(*imagecv, weightcv);
  // should be around 1....
  cout << "normalizedSig of convolved image " <<  cvImageNormalizedSig << endl;
  // write it somewhere for later use
  SaveNormalizedSig();

  
  // normalize imagecv : imagecv *= sqrt(weightcv);
  weightcv.ApplyFun(sqrt);
  *imagecv *= weightcv;

  // write image for studies 
  if (getenv ("SAVECV")) { FitsImage toto("savecv.fits",*imagecv);}
  
  Pixel mean, sig;
  imagecv->SkyLevel(imFrame, &mean, &sig);
  cvImageNormalizedMean = mean;
  cout << " compute the mean and sigma of img*sqrt(weight) : "  
       << mean << ' ' << sig << endl;
  
}



#define NORMALIZED_SIG_KEY "NORMSIG"

void DetectionProcess::SaveNormalizedSig() const
{
  FitsHeader headw(weightName, RW);
  headw.AddOrModKey(NORMALIZED_SIG_KEY,cvImageNormalizedSig," average error excess of PSF photometry ");
}

bool DetectionProcess::GetNormalizedSig()
{
  FitsHeader headw(weightName);
  if (headw.HasKey(NORMALIZED_SIG_KEY))
    {
      cvImageNormalizedSig = headw.KeyVal(NORMALIZED_SIG_KEY);
      return true;
    }
  else return false;
}

    

void DetectionProcess::FirstDetections(DetectionList &List) const
{
  Pixel threshold = cvImageNormalizedMean + nSigDet *  cvImageNormalizedSig;
  int radLocmax = int(radAper*0.5)+1;
  for (int j=0; j<image->Ny(); ++j)
    for (int i=0; i<image->Nx(); ++i)
      {
	if ((*imagecv)(i,j) > threshold 
	    && IsLocalMax(*imagecv, i, j, radLocmax))
	  {
	    List.push_back(new Detection(double(i), double(j),(*imagecv)(i,j) ));
	    // save the original sig2Noise
	    Detection *det = List.back();
	    det->sig2NoiseCv = (*imagecv)(i,j); 
	  }
      }
}
	  


void DetectionProcess::DetectionPosition(Detection &Det) const
{
  int tx = image->Nx() ;
  int ty = image->Ny() ;
  int d = (int) (floor(radWeight+1));
  double rad_weight2 = sqr(radWeight) ;
  double rad_bad2 = sqr(radBad);
  double sx=0;
  double sy=0;
  double s=0;
  double SumWi = 0, Wi = 0;
  int ix = (int) (floor(Det.x)) ;
  int iy = (int) (floor(Det.y)) ;
  double back = LocalImageBackground(ix, iy);
  Det.LocalBack() = back ;

  double dist_bad_init = sqr(1000);
  Det.distBad = dist_bad_init;
  for (int l=-d; l < d+1; l++)
    {
      for (int m=-d; m < d +1;  m++)
	{
	  int xu = ix + l ; 
	  int yu = iy + m ;
	  double r2=sqr(double(xu) - Det.x) + sqr(double(yu)-Det.y);
	  if ( ( xu < tx ) &&  (xu >= 0 ) &&
	       ( yu < ty ) &&  (yu >= 0 ) && (r2 <rad_weight2) )
	    {
	      double val = (*image)(xu,yu)- back ;
	      double localWeight = (*weight)(xu,yu);
	      if (localWeight == 0)
		{
		  if ( r2 < rad_bad2) Det.nBad++;
		  Det.distBad = min(Det.distBad, r2);
		}
	      Wi = exp(-(r2)/(2*sigFilter*sigFilter))*localWeight;
	      s += val*Wi ;
	      sx += l*val*Wi;
	      sy += m*val*Wi;
	      SumWi += Wi;
	    }
	}
    }
  if (fabs(SumWi) > 1.e-20)
    Det.flux = s/SumWi ; // very bizarre
  Det.distBad = (Det.distBad == dist_bad_init)? -1 : sqrt(Det.distBad);
  if (s != 0) 
    {
      sx  = sx/s;
      sy  = sy/s;
      Det.x = sx + ix;
      Det.y = sy + iy;
     }
  return; 
}


void DetectionProcess::RefineDetectionPosition(Detection &Det) const
{
  int count=0;
  double deltax, deltay;
  // xs and ys are the "starting position". save them for posterity
  Det.xs = Det.x;
  Det.ys = Det.y;
  do 
    {
      double xbar = Det.x;
      double ybar = Det.y;
      DetectionPosition(Det);
      deltax = fabs(Det.x-xbar);
      deltay = fabs(Det.y-ybar);
      count++;
      if (count>10) break;
    } while (deltax>0.005 || deltay>0.005);
}

// computes everything but the position
void DetectionProcess::SetDetectionScores(Detection &Det) const
{
#ifdef DEBUG
  //cout << "DetectionProcess::SetDetectionScores" << endl;
#endif
  int tx = image->Nx();
  int ty = image->Ny();
  int d = (int) ((radWeight + 1.));
  double rad_weight2 = sqr(radWeight);
  double rad_aper2 = sqr(radAper);
  // accumulators
  double sx=0;
  double sy=0;
  double S=0;
  double SumWi = 0, SumSigx = 0, SumSigy = 0, SumSigxy = 0;
  int ix = (int) (Det.x + 0.5 ) ; // nearest integer
  int iy = (int) (Det.y + 0.5 ) ;
  double back = LocalImageBackground(ix, iy);
  Det.LocalBack() = back ;

  double aperFlux = 0; double aperFluxDeno = 0;
  /* recompute flux in this routine, because we want to use it for 
     measuring fluxes ob object detected by other means. */
  double flux = 0; double fluxWeight = 0;
  double varXX = 0, varYY = 0, varXY = 0;
  double sumPsf = 0;
  int area = 0;
  int bad = 0;
  for (int l=-d; l <= d; l++)
    {
      for (int m=-d; m <=d ;  m++)
        {
          int xu = ix + l ;
          int yu = iy + m ;
	  double dx = xu - Det.x;
	  double dy = yu - Det.y;
          double dist2 = sqr(dx) + sqr(dy);
          if ( ( xu < tx ) &&  (xu >= 0 )&&
               ( yu < ty ) &&  (yu >= 0 ))
            { 
	      const double localWeight = (*weight)(xu, yu);
	      const double val = (*image)(xu,yu)- back;
	      if (dist2 <rad_weight2) 
		{

		  double psfVal  = exp(-(dist2)/(2*sigFilter*sigFilter));
		  sumPsf += psfVal; // sum psf
		  double Wi = psfVal*localWeight;
		  // flux
		  flux += Wi*val; // sum (psf * w * val)
		  fluxWeight += psfVal * Wi; // sum psf^2*w
		  // shape params
		  sx += dx*val*Wi;
		  sy += dy*val*Wi;
		  SumSigx += sqr(dx) * val * Wi;
		  SumSigy += sqr(dy) * val * Wi;
		  SumSigxy += dx * dy * val* Wi;
		  SumWi += Wi;
		  S += Wi * val;
		  // position variance
		  varXX += sqr(dx) * Wi * psfVal;
		  varYY += sqr(dy) * Wi * psfVal;
		  varXY += dx * dy * Wi * psfVal;
		  // varPosDeno += Wi * val; already compudted in flux
		}
	      if (dist2 < rad_aper2)
		{
		  aperFlux += val*localWeight;
		  aperFluxDeno += localWeight;
		  area++;
		  if(localWeight == 0.0) bad++;
		}
            }
        }
    }
  if (fluxWeight > 0)
    {
      Det.flux = flux*sumPsf/fluxWeight;
      //      Det.eFlux^2 = sqr(sumPsf)/fluxWeight;
      Det.eFlux = sumPsf/sqrt(fluxWeight);
      // The above formula assumes that the pixels are uncorrelated,
      // which is wrong. Correct by the average factor:
      Det.eFlux *= cvImageNormalizedSig;
      Det.sig2Noise = Det.flux/Det.eFlux;
    }
  else Det.eFlux = -1;
  if (flux != 0)
    { // no need to account for sumPSf != 1
      Det.vx = varXX/sqr(flux);
      Det.vy = varYY/sqr(flux);
      Det.vxy = varXY/sqr(flux);
    }
  if ( fabs(S) > 1.e-20 )
    {
      Det.Mxx() = SumSigx / S - sqr(sx/S);
      Det.Myy() = SumSigy / S - sqr(sy/S);
      Det.Mxy() = SumSigxy / S - sx*sy/sqr(S);
    }
  Det.aperFlux = (aperFluxDeno >0)? area*aperFlux/aperFluxDeno : 0;
  Det.Area() = area;
  Det.NBad() = bad;
}

#include "vutils.h" // for FarrayMedian

double DetectionProcess::LocalImageBackground(const int i, const int j) const
{
   if (radBackMax == 0) return 0;
  int tx = image->Nx();
  int ty = image->Ny();

  static Pixel *pixels = NULL;
  static int stampSize = 0; // the size with which it was reserved
  int backSide = int (radBackMax+1.); // pour etre sur que backSide  > radBackMin au cas ou l'anneau est fin ...
  if (pixels == NULL || backSide > stampSize)
    {
      if (pixels) delete [] pixels;
      pixels = new Pixel[(2*backSide+1)*(2*backSide+1)];
      stampSize = backSide;
    }
  int pixelCount = 0;
  double distCut2 = sqr(radBackMin);
  for (int l=-backSide; l <= backSide; l++)
    {
      for (int m=-backSide; m <=backSide ;  m++)
        {
          double dist2 = sqr(m)+sqr(l);
	  // forget pixels in the inner circle
	  if (dist2 < distCut2) continue;
          int xu = i + l ;
          int yu = j + m ;
          if ( ( xu < tx ) &&  (xu >= 0 )&&
               ( yu < ty ) &&  (yu >= 0 ))
            { 
	      double localWeight = (*weight)(xu, yu);
	      if (localWeight == 0) continue;
	      pixels[pixelCount++] = (*image)(xu,yu);
	    }
	}
    }
  return FArrayMedian(pixels, pixelCount); 
}


static bool IgnoreIt(const BaseStar *star )
{
  return (star->flux == 0);
}

#include "fastfinder.h"
/* The new cleaner routine uses the fastfinder.  It looks for the
nearest star from the current star, if the neighbour is in the cut
radius and have a flux greater than the currrent star, then the
current star is set to zero and then deleted from the list. This trick
is needed because the fastfinder cannot accomodate stars which
disappear.*/

static void CleanDetectionList(DetectionList & stl, double delta)
{
  cout << "Size of the list before cleaning: " << stl.size() << endl ;
  cout << "Cut radius for double detectiom " << delta << endl;
  // sort and copy of the list for old cleaner
  stl.FluxSort() ;
  
  FastFinder finder(*((BaseStarList*) &stl));
  
  for (DetectionIterator it = stl.begin(); it != stl.end(); it++)
    {
      Detection *star = *it ;
      // To avoid having neighbour = current star using FindClosest
      // The flux is reaffected to its initial value after.
      double KeepFlux = star->flux;
      star->flux = 0;
      bool keep_it = true;
      const BaseStar *neighbour = finder.FindClosest(*star, delta, &IgnoreIt); 
      if ( neighbour != NULL)
	{
	  if ( star->flux < neighbour->flux)
	    keep_it = false;
	}
      if (keep_it) star->flux = KeepFlux;
    }
  
  // to erase the star with 0 flux.
  for (DetectionIterator it = stl.begin(); it != stl.end();)
    {
      Detection * star = *it ; 
      if ( star->flux == 0)
	it = stl.erase(it);
      else ++it;
    }
  return ;
}


void DetectionProcess::FluxFromPos(const Point &Where, double &Flux)
{
  Detection det(Where.x, Where.y);
  if (!weight) weight = new FitsImage(weightName);
  SetDetectionScores(det);
  Flux = det.flux;
}
  

void DetectionProcess::DoDetection(DetectionList &List)
{
  cout << " starting detection on " << imageName << endl;
  SetCuts(PrintDatacards + PrintCuts);
  ImagesConvolve();
  FirstDetections(List);
  // do not nee any longer imagecv
  delete imagecv; imagecv = NULL;

  cout << " number of raw detections " << List.size() << endl;


  //need weight
  weight = new FitsImage(weightName);
  // refine detections, compute flux, and number of bad pixels.
  Frame imageBoundaries(*weight); // this is the total area of the image
  for (DetectionIterator i = List.begin(); i != List.end(); )
    {
      Detection &det = *(*i);
      RefineDetectionPosition(det);
      // to remove objects that drifted outside the image
      if (!imageBoundaries.InFrame(det))
	i = List.erase(i);
      else ++i;
    }
  cout << " number after position refinement and removal of drifters " 
       << List.size() << endl;

  // clean stuff that contains pixels with null weight
  for (DetectionIterator i = List.begin(); i != List.end(); )
    {
      Detection &det = **i;
      if (det.NBad()) i = List.erase(i);
      else ++i;
    }

  cout << " number of detections after removal of bad " << List.size() << endl;

  // remove identical stuff.
  CleanDetectionList(List, minDistDouble);

  cout << " number of detections after removal of neighbours " 
       << List.size() << endl;
  
  // fill the scores
  for (DetectionIterator i = List.begin(); i != List.end(); ++i)
    SetDetectionScores(*(*i));

  for (DetectionIterator i = List.begin(); i != List.end(); )
    {
      Detection &det = **i;
      if (det.sig2Noise < nSig)
	i = List.erase(i);
      else ++i;
    }

  cout << " number of detections after cut on S/N : " 
       << List.size() << endl;

}


void DetectionProcess::DetectionScoresFromPositions(const BaseStarList & Pos, 
						  DetectionList &Out,
						  const bool FixedPositions)
{
  /* we need a value for cvImageNormalizedSig
     if the image alread had a detection  run on it it is in the header,
     if not, we have to compute it */
  cout << " starting addressed detection on " << imageName << endl;
  cout << "  number of given positions :" << Pos.size() << endl;
  SetCuts(PrintCuts);
  if (!GetNormalizedSig())
    {
      cout << " have to convolve image to get the actual noise " << endl;
      ImagesConvolve();
      if (imagecv) delete imagecv; imagecv = NULL;
    }
  else
    {
      cout << " Normalized sigma of convolved images " 
	   << cvImageNormalizedSig << endl;
    }

  // weight is not loaded by constructor.
  weight = new FitsImage(weightName);

  Out.clear();
  for (BaseStarCIterator in = Pos.begin(); in != Pos.end(); ++in)
    {
      const BaseStar &inPos = **in;
      Detection *out = new Detection(inPos.x, inPos.y);
      if (!FixedPositions) RefineDetectionPosition(*out);
      SetDetectionScores(*out);
      Out.push_back(out);
    }
  if (Pos.size() != Out.size())
    {
      cerr << " ERROR : Serious bug in DetectionScoresFromPositions " << endl
	   << " input and output list have different sizes " 
	   << Pos.size() << ' ' << Out.size() << endl;
    }
}
      

#include "reducedimage.h"
#include "wcsutils.h"
#include "gtransfo.h"
#include "sestar.h"

void DetectionProcess::SetScoresFromRef(DetectionList &List, 
					const ReducedImage &Ref)
{
  SetCuts(PrintCuts);
  SEStarList RefList(Ref.CatalogName());
  FastFinder finder(*SE2Base(&RefList));
  Gtransfo *Pix2RaDec = NULL;
  WCSFromHeader(imageName, Pix2RaDec);
  for (DetectionIterator i= List.begin(); i!= List.end(); ++i)
    {
      Detection &Det = **i;
      if (Pix2RaDec) // some simulations don't have WCS's. no else needed
	{
	  FatPoint raDec;
	  Pix2RaDec->TransformPosAndErrors(Det, raDec);
	  Det.vRaRa = raDec.vx;
	  Det.vDecDec = raDec.vy;
	  Det.vRaDec = raDec.vxy;
	}
      FluxFromPos(Det, Det.fluxRef);
      /* put a positive prctInc even when fluxRef<0. this simplifies 
	 quality cuts... */
      if (Det.fluxRef>0) Det.prctInc = Det.flux/Det.fluxRef;
      else Det.prctInc = 100;
      const SEStar *neighbour = (const SEStar *) finder.FindClosest(Det,20.);
      if (neighbour)
	{
	  Det.fluxObjRef = neighbour->flux;
	  Det.xObj = neighbour->x;
	  Det.yObj = neighbour->y;
	  Det.distObjRef = Det.Distance(*neighbour);
	}
    }
}



/************************ Detection *******************************/


Detection::Detection(const double X, const double Y, const double Flux) 
  : BaseStar(X,Y,Flux)
{
  eFlux = 0;
  sig2Noise = 0;
  aperFlux = 0;
  sig2NoiseCv = 0;
  xs = 0;
  ys = 0;
  vx = 0; 
  vy = 0; 
  vxy = 0; 
  area = 0; 
  nBad = 0; 
  distBad = 0; 
  mxx = 0; 
  myy = 0;
  mxy = 0;
  //
  ra = dec = 0;
  vRaRa =  vDecDec = vRaDec = 0;
  fluxRef = 0;
  prctInc = 0;
  xObj = yObj = 0;
  fluxObjRef = distObjRef = 0;
  localback = 0. ;
}


#include "fastifstream.h"

void Detection::read_it(fastifstream& r, const char * Format)
{
  int format = DecodeFormat(Format, "Detection");
  BaseStar::read_it(r, Format);
  r >> eFlux 
    >> sig2Noise 
    >> aperFlux 
    >> sig2NoiseCv 
    >> xs >> ys;
  if (format<3) // moved position uncertainties into BaseStar 
    r >> vx >> vy >> vxy;
  r >> area 
    >> nBad 
    >> distBad 
    >> mxx >> myy >> mxy
    >> ra >> dec 
    >> vRaRa >> vDecDec >> vRaDec
    >> fluxRef
    >> prctInc
    >> xObj >> yObj
     >> fluxObjRef >> distObjRef
    ;
  /* the test that comes just after is far from being full proof: the
     actual format was changed w/o changing the label. What shows the
     format change is the presence of "fwhmref" in the file header. 
     We don't have a routine that checks for that yet, and XML I/Os 
     make this question obsolete.
  */
  if (format < 2)
    {
      double fwhmRef, shapeRef;
      r >> fwhmRef >> shapeRef;
    }
  if (format >=1)
    r >> localback ;
}


Detection* Detection::read(fastifstream& r, const char* Format)
{
  //awful stuff (we do not have a persistence system..)
  if (strstr(Format,"MatchedDetection"))
    {// degrade it to a Detection
      MatchedDetection *md = MatchedDetection::read(r,Format);;
      Detection *d = new Detection;
      *d = *md;
      delete md;
      return d;
    }
  else
    {
      Detection *d = new Detection;
      d->read_it(r, Format);
      return d;
    }
}

std::string Detection::WriteHeader_(ostream & stream, const char*i) const
{
  if (i==NULL) i = "";
  string baseStarFormat =  BaseStar::WriteHeader_(stream, i);
  stream << "# eflux"<< i <<" : flux error" << endl
	 << "# sig2noi"<< i <<" : " << endl
  	 << "# apflux"<< i <<" : aperture flux" << endl
  	 << "# s2ncv"<< i <<" : sig2noise at detection" << endl
	 << "# xs" << i << " : x detection position " << endl
	 << "# ys" << i << " : y detection position " << endl
    /* position variance now in BaseStar (FatPoint indeed)
  	 << "# varxx"<< i <<" : pos variance" << endl
  	 << "# varyy"<< i <<" : pos variance" << endl
  	 << "# varxy"<< i <<" : pos variance" << endl
    */
  	 << "# area"<< i <<" : number of pixels used in aperflux" << endl
  	 << "# nbad"<< i <<" : number of bad pixels use in flux" << endl
	 << "# distbad" << i << " : distance to nearest bad pix" << endl
    	 << "# mxx"<< i <<" : shape param" << endl 
    	 << "# myy"<< i <<" : shape param" << endl 
    	 << "# mxy"<< i <<" : shape param" << endl
	 << "# ra"<< i <<" : RA" << endl
	 << "# dec"<< i <<" : dec" << endl
	 << "# vrara" << i<<" : var(ra)" << endl
	 << "# vdede" << i<<" : var(dec)" << endl
	 << "# vrade" << i<<" : covar(va,dec)" << endl
	 << "# fluxre" << i<<" : flux on the ref under SN" << endl
	 << "# prctinc" << i<<" : flux/fluxref if fluxref>0, 100 if not" << endl
	 << "# xobj"<<i<<" : x pos of nearest object " << endl
	 << "# yobj"<<i<<" : y pos of nearest object " << endl
	 << "# fobjref"<<i<<" : flux of nearest object " << endl
	 << "# distobj"<<i<<" : distance to nearest object " << endl
    //	 << "# fwhmref" <<i << " : FWHM of the nearest object"  << endl
    // << "# shaperef" << i << " : -2.5*log10(flux/fluxmax) of the nearest objet" << endl
	 << "# back"<<i<<" : local background " << endl
	      
    ;
  /*
    format == 3 means that varXX, varYY, and varXY just disappear : 
    we now use the ones in BAseStar (FatPoint indeed) */
  return baseStarFormat + " Detection 3 "; 
}

#include <iomanip>

void Detection::writen(ostream &s) const 
{
  ios::fmtflags  old_flags =  s.flags();
  BaseStar::writen(s);
  s  << eFlux << ' ' //27
    << sig2Noise << ' ' 
    << aperFlux << ' ' 
    << sig2NoiseCv << ' ' 
    << xs << ' ' << ys << ' '
    /* 
       << varXX << ' ' << varYY << ' ' << varXY << ' ' 
       position uncertainties moved into BaseStar 
    */
    << area << ' ' 
    << nBad << ' ' 
    << distBad << ' ' 
     << mxx << ' ' << myy << ' ' << mxy << ' '; // 13 14 15
  int prec = s.precision();
  s << setprecision(10);
  s << ra << ' ' << dec << ' ';
  s << setprecision(prec);
  s << vRaRa << ' ' << vDecDec << ' '  << ' ' << vRaDec << ' '
    << fluxRef << ' '
    << prctInc << ' '
    << xObj << ' ' << yObj << ' ' 
    << fluxObjRef << ' ' << distObjRef
    << ' ' << localback << ' '
    ;
  s.flags(old_flags);
}

  

#include "starlist.cc"
template class StarList<Detection>; // to instantiate

BaseStarList *Detection2Base(DetectionList *D)
{
  return (BaseStarList *) D;
}


/************* MatchedDetection ***************/

std::string MatchedDetection::WriteHeader_(ostream & stream, 
					   const char*i) const
{
  
  if (i==NULL) i = "";
  string format = Detection::WriteHeader_(stream,i);
  stream << "# count"<< i <<" : number of associated detections " << endl;
  for (unsigned k=1; k<=others.size() ; ++k)
    {
      char st[8];
      sprintf(st,"%s_%d",i,k);
      format += Detection::WriteHeader_(stream, st);
    }
  return " MatchedDetection 0 " + format;
}


MatchedDetection* MatchedDetection::read(fastifstream& r, const char* Format)
{
  const char *subFormat = strstr(Format," Detection");
  Detection *d = Detection::read(r, subFormat);
  MatchedDetection *result = new MatchedDetection(*d);
  delete d;
  unsigned count;
  r >> count;
  for (unsigned k=0; k<count ; ++k)
    {
      result->others.push_back(Detection::read(r, subFormat));
    }
  return result;
}

void MatchedDetection::writen(ostream &s) const 
{
  Detection::writen(s);
  s << others.size() << ' ';
  for (unsigned k=0; k<others.size(); ++k)
    {
      others[k]->writen(s);
    }
}

bool MatchedDetection::CompatibleSig2Noise(const double &Fact) const
{
  if (Fact < 0 ) return true ;
  for (unsigned k=0; k<others.size(); ++k)
    {
      if (Fact*others[k]->Sig2Noise() < Sig2Noise())
	  return false;
    }
  return true;
}

bool MatchedDetection::CompatiblePosition(const double &DistMax) const
{
  for (size_t k=0; k<others.size(); ++k)
    if (Distance(*others[k]) < DistMax)
      return false;
  return true;
}

/************** MatchedDetectionList *************/

MatchedDetectionList::MatchedDetectionList(const DetectionList &L)
{
  for (DetectionCIterator i = L.begin(); i!=L.end(); ++i)
    push_back(new MatchedDetection(**i));
}

bool MatchedDetectionList::OneToOneAssoc(const string &ImageName, 
					 DetectionList &L)
{
  if (size() != L.size())
    {
      cerr << " Cannot associate 1 to 1 detections from " << ImageName
	   << " because existing list has a different size " << endl;
      return false;
    }
  DetectionIterator di = L.begin();
  for (MatchedDetectionIterator mi = begin(); mi != end(); ++mi, di++)
    {
      (*mi)->AssociateDet(*di);
    }
  imageNames.push_back(ImageName);
  return true;
}


void MatchedDetectionList::ApplyCuts()
{
  DatDet datdet(DefaultDatacards());
  cout << "On MatchedDetectionList (" << size() << ") :" <<  endl;
  cout << " appliying compatibility cut on partial subs : "
       <<" factor =" << datdet.compFactor << endl;
  for (MatchedDetectionIterator i = begin(); i != end(); )
    {
      MatchedDetection &md = **i;
	if (!md.CompatibleSig2Noise(datdet.compFactor))
	  i = erase(i);
	else ++i;
    }
  
  cout << "New size of MatchedDetectionList : " << size() <<  endl;
}  


#include "starlist.cc"
template class StarList<MatchedDetection>; // to instantiate

