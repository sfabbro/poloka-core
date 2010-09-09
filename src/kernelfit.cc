#include <math.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <assert.h>
#include <fstream>
#include <memory> // auto_ptr

#include "kernelfit.h"
#include "basestar.h"
#include "vutils.h"
#include "matvect.h"
#include "starmatch.h"

#include "iohelpers.h"
#include "polokaexception.h"
#include "scorescollection.h"
#include "imagepair.h"
#include "reducedutils.h" // for MedianPhotomRatio
#include "fastfinder.h"

/*
  ToDo :

  - add contribution of object to weight (AssignWeight) (using Gaussian model ?)
        - Need gains to do that..... 
  - handle stamps which are not fully inside images (in extract_pixels)
        idea :  shift center accordingly


*/


static double sqr(double x) { return x*x;}

/**************    Stamp and StampList **********************/

//! an odd size DImage extracted from an Image, centered on (xc,yc)
class Stamp {  
  /*
private:
  int hBestSize;
  int hWorstSize;
  */
public :

  //! center in the source Image
  int xc,yc ;  

  //! pixels (in double precision until one changes DPixel definition )
  DImage bestPixels;

  //! cookup for weights used when fitting the kernel
  DImage weight;

  //!
  DImage bestVariance,worstVariance;

  //!
  double cutBest, cutWorst;

  int nActivePix;

  //! will be filled and used to discard outliers from the fit 
  double chi2; 
  //  int HSize () const {return hsize;}  

  //! we should not assume that there is a BaseStar corresponding to this stamp, but if any, put its pointer here
  BaseStarRef star; 
  Stamp(const BaseStar *Star, const Image& image, 
	const int HBestSize, const int HWorstSize);

  int AssignWeight(ImagePair &ImPair, const Kernel &GuessedKernel, const StarMatch &Sm);


  int UpdateWeight(const Kernel &K);


};

//! Nothing but a list of Stamps
class StampList : public list<Stamp>
{
public :
  StampList(ImagePair &ImPair, const StarMatchList &List, 
	    const int hStampSize, const Kernel &GuessedKernel, const int MaxStamps);

  void init(ImagePair &ImPair, const StarMatchList &List, 
	    const int hStampSize, const Kernel &GuessedKernel, const int MaxStamps);


  double Sig2Noise() const;

};

typedef list<Stamp>::iterator StampIterator;
typedef list<Stamp>::const_iterator StampCIterator;

//***********************  Stamp ***************************

static void extract_pixels(DImage &Target, const Image &image, 
			   const int xc, const int yc)
{
  int hStampSize = Target.Nx()/2;
  int xoff = xc - hStampSize;
  int yoff = yc - hStampSize;
  int xstart = max(0,xc - hStampSize);
  int ystart = max(0,yc - hStampSize);
  int width = Target.Nx();
  int height = Target.Ny();
  int xend = min(xstart+width,image.Nx());
  int yend = min(ystart+height, image.Ny());
  if (width*height != (xstart-xend)*(ystart-yend))
    {
      cerr << " we miss pixels for (" << xc << "," << yc << ")" << endl;
      abort(); // added it in feb 2010, to see if it ever shows up
    }
  for (int j=ystart; j <yend; ++j)
    {
      for (int i= xstart; i< xend; i++)
	{
	  Target(i-xoff,j-yoff) = image(i,j);
	}
    }
}

#include "apersestar.h"


Stamp::Stamp(const BaseStar *Star, const Image& image, const int HBestSize, const int HWorstSize)
  : xc(int(round(Star->x))), yc(int(round(Star->y))), nActivePix(0), star(Star)
{
  bestPixels.Allocate(2*HBestSize+1,2*HBestSize+1);
  weight.Allocate(2*HWorstSize+1, 2*HWorstSize+1);
  extract_pixels(bestPixels, image, xc,yc);
}


static void add_object_noise(const BaseStar *B, 
			     const int xc, const int yc,
			     const double FluxOverGain,
			     DImage &Variance)
{
 double mxx, myy, mxy;
  const AperSEStar* ap = dynamic_cast<const AperSEStar*>(B);
  if (ap)
    {
      mxx = ap->gmxx;
      myy = ap->gmyy;
      mxy = ap->gmxy;
    }
  else 
    {
      const SEStar* se = dynamic_cast<const SEStar*>(B);
      assert(se != NULL);
      mxx = se->Mxx();
      myy = se->Myy();
      mxy = se->Mxy();
    }
  // account for object noise using an object "model".
  int hsize = Variance.Nx()/2;
  double dx = B->x-xc+hsize;
  double dy = B->y-yc+hsize;
 
  double det = mxx*myy-sqr(mxy);
  double wxx = myy/det;
  double wyy = mxx/det;
  double wxy = -mxy/det;
  double factor = FluxOverGain/(2*M_PI*sqrt(det));

  for (int j=0; j< Variance.Ny(); ++j)
    {
      double y = j-dy;
      for (int i=0; i< Variance.Nx(); ++i)
	{
	  double x = i-dx;
	  Variance(i,j) += factor*exp(-0.5*(wxx*x*x+wyy*y*y+2*wxy*x*y));
	}
    }
}
  
void weight_to_variance(DImage &Weight, int &NPix, double &VAverage)
{
  // compute the average weight 
  double waverage = 0;
  NPix = 0;
  const DPixel *pend = Weight.end();
  for (DPixel* p = Weight.begin(); p < pend; ++p)
    if (*p) { waverage += (*p); NPix++;}  
  if (!NPix) return;
  VAverage = NPix/waverage; // average variance
  // transform to variance
  for (DPixel* p = Weight.begin(); p < pend; ++p) 
    if (*p) (*p) = 1./(*p); else (*p) = 1E15*VAverage;
}


static int build_variance(const Image &ImageWeight, 
			  const double &Xc, const double &Yc,
			  const double &Gain, const BaseStar *Star, 
			  DImage &VarStamp, double &CutVar)
{
  extract_pixels(VarStamp, ImageWeight, Xc, Yc);
  // VarStamp constains weights
  int count;
  double vAverage;
  weight_to_variance(VarStamp, count, vAverage);
  if (!count) return 0; // we are done: the stamp is dead
  CutVar = 1e8*vAverage;

   // tmpweight is now a variance
  add_object_noise(Star, Xc, Yc, Star->flux/Gain, VarStamp);
  return count;
}
  



//! returns the number of non zero pixels
int Stamp::AssignWeight(ImagePair &ImPair,  
			const Kernel &GuessedKernel, 
			const StarMatch &Sm)
{
  const Image &bestWeight = ImPair.BestWeight();
  const Image &worstWeight = ImPair.WorstWeight();
  double bestGain = ImPair.BestGain();
  double worstGain = ImPair.WorstGain();

  bestVariance.Allocate(bestPixels.Nx(),bestPixels.Ny());
  int cBest = build_variance(bestWeight, xc, yc, bestGain, Sm.s1, 
			    bestVariance, cutBest);
  if (cBest == 0) return 0;

  worstVariance.Allocate(weight.Nx(),weight.Ny());
  int cWorst = build_variance(worstWeight, xc, yc, worstGain, Sm.s2,
			     worstVariance, cutWorst);

  if (cWorst ==0) return 0;
  return UpdateWeight(GuessedKernel);

}

// compute Var(K*B+W) = k2*Var(B) + Var(W) and inverts it.
int Stamp::UpdateWeight(const Kernel &K)
{
  Kernel K2(K);
  K2 *= K; // reminder: Var(a*x) = a**2 Var(x)
  // weight and  bestVariance are variances. weight is allocated in Stamp constructor.
  Convolve(weight, bestVariance, K2); 

  nActivePix = 0;
  for (int j=0; j<weight.Ny();  ++j)
    for (int i=0; i<weight.Nx(); ++i)
    {
      double &val = weight(i,j);// contains bestVariance (convolved)
      double worstVar = worstVariance(i,j);
      if (val > cutBest || worstVar > cutWorst)  val =0;
      else
	{
	  val = 1./(val + worstVar);
	  nActivePix++;
	}    
    }

  // zero infinite variances (for noise computation)
  for (int j=0; j<bestVariance.Ny();  ++j)
    for (int i=0; i<bestVariance.Nx(); ++i)
      {
	double &val = bestVariance(i,j);
	if (val > cutBest) val = 0;
      }

  return nActivePix;
}


StampList::StampList(ImagePair &ImPair, const StarMatchList &Objects, 
		     const int hBestSize, const Kernel &GuessedKernel, const int MaxStamps)
{
  init(ImPair, Objects, hBestSize, GuessedKernel, MaxStamps);
}

void StampList::init(ImagePair &ImPair, const StarMatchList &Objects, 
		     const int hBestSize, const Kernel &GuessedKernel, const int MaxStamps)
{
  int count = 0;
  const Image& bestImage  = ImPair.BestImage();
  // What are the image boundaries ? This is  the hard limit :
  const Frame& imageFrame(bestImage);
  // we could consider the more stringent the soft limit :
  //   const Frame& imageFrame = ImPair.CommonFrame();
  int hWorstSize = hBestSize-GuessedKernel.HSizeX();

  int nPixCut = (2*hWorstSize+1)*(2*hWorstSize+1)/2;

  for (StarMatchCIterator sm = Objects.begin(); sm != Objects.end() && count < MaxStamps; ++sm)
    {
      const BaseStar *s = sm->s1;
      if (imageFrame.MinDistToEdges(*s) < hBestSize+1) continue;
      Stamp stamp(s, bestImage, hBestSize, hWorstSize);
      // compute a plausible weight image for this stamp
      if (stamp.AssignWeight(ImPair, GuessedKernel, *sm) > nPixCut )
       {
	 push_back(stamp);
	 ++count;
       }
    }
}

double StampList::Sig2Noise() const
{
  double signal2 = 0;
  double noise2 = 0;
  for (StampCIterator i = begin(); i != end(); ++i)
    {
      const BaseStar *b = i->star;
      const SEStar* se = dynamic_cast<const SEStar*>(b);
      signal2 += sqr(se->flux);
      noise2 += sqr(se->EFlux());
    }
  return (noise2>0) ? sqrt(signal2/noise2) : 0;
}
      


/**************  end of  Stamp and StampList **********************/

/****************************   selection of the objects for fit ***************************************/
static bool DecreasingFluxMax(const SEStar *S1, const SEStar *S2)
{
  return (S1->Fluxmax() > S2->Fluxmax());
}

// another starlist selection routine
static bool GoodForFit(const SEStar *Star, const double &SaturLev, const double &Bmin, const double& minsignaltonoiseratio)
{
  return ((Star->FlagBad() == 0 ) && 
	  ( Star->Flag() <= 3 )   && //(keep blended objects) 
	  ( Star->Fluxmax()+Star->Fond() < SaturLev ) &&
	  ( Star->B() > Bmin ) &&
	  ( Star->Flux_auto()/Star->Eflux_auto() > minsignaltonoiseratio)
	  );
}

static int MakeObjectList(const ImagePair& ImPair, const double MaxDist,
			  StarMatchList &objectsUsedToFit)
{
  //get list of objects suitable for the kernel fit 
  const ReducedImageRef &best= ImPair.Best();
  const ReducedImageRef &worst= ImPair.Worst();
  const Frame &intersection = ImPair.CommonFrame();


  AperSEStarList bestStarList(best->AperCatalogName());
  
  // check whether background was subtracted
  // TODO : check if this is useful by any mean.
  if(best->BackSub()) {
    // then set all fond() to zero
    SetStarsBackground((SEStarList &)bestStarList,0.);
  }
  
  // sort with decreasing quality for the kernel fit. The chosen criterion is the peak flux (e.g. in best)
  bestStarList.sort(DecreasingFluxMax);


  cout << " Entering FilterStarList with "<< bestStarList.size() 
       << " stars " << endl;
  //star list selection
  SEStarList worstStarList(worst->ImageCatalogName());
  
  // check whether background was subtracted
  if(worst->BackSub()) {
    // then set all fond() to zero
    SetStarsBackground(worstStarList,0.);
  }

  double satfactor = 0.95;
  double saturLevBest = best->Saturation() * satfactor;
  double saturLevWorst = worst->Saturation() * satfactor;  
  double bfactor = 0.2;
  double bMinBest = best->Seeing() * bfactor;
  double bMinWorst = worst->Seeing()* bfactor;
  double mindist = worst->Seeing()*3;
  double minsignaltonoiseratio = 10;
  cout << "cuts for best, satur bmin " << saturLevBest << " " << bMinBest << endl;
  cout << "cuts for worst, satur bmin " << saturLevWorst << " " << bMinWorst << endl;
  FastFinder worstFinder(*SE2Base(&worstStarList));

  int n_good_for_fit_best  = 0;
  int n_matches = 0;
  int n_good_for_fit_worst = 0;
  int n_in_frame = 0;
  int n_far_from_edges = 0;

    
  for (AperSEStarIterator sibest = bestStarList.begin(); sibest != bestStarList.end(); sibest++)
    {
      AperSEStar *sb = *sibest;
      const SEStar *sworst;
      if (sb->gflag == 0 && GoodForFit(sb, saturLevBest, bMinBest, minsignaltonoiseratio)) {
	n_good_for_fit_best++;
	sworst = (SEStar *) worstFinder.FindClosest(*sb,MaxDist);
	if (sworst) {
	  if(sb->Distance(*sworst) < MaxDist) { // useless?
	    n_matches++;
	    if(GoodForFit(sworst, saturLevWorst, bMinWorst, minsignaltonoiseratio)) {
	      n_good_for_fit_worst++;
	      if(intersection.InFrame(*sb)) {
		n_in_frame++;
		if(intersection.MinDistToEdges(*sb) > mindist) // remove objects close to the edge
		  {
		    n_far_from_edges++;
		    objectsUsedToFit.push_back(StarMatch(*sb,*sworst,sb,sworst));
		  }
	      }
	    }
	  }
	}
      }
    }
  // dump all counters
  cout << "n_good_for_fit_best " << n_good_for_fit_best << endl;
  cout << "n_matches " << n_matches << endl;
  cout << "n_good_for_fit_worst " << n_good_for_fit_worst << endl; 
  cout << "n_in_frame " << n_in_frame << endl;
  cout << "n_far_from_edges " << n_far_from_edges << endl; 
  
  return objectsUsedToFit.size();
}


/***********************   End of object list making **************************/




//#define DEBUG

#define DO2(A) A;A;
#define DO5(A) A;A;A;A;A;
#define DO10(A) DO2(DO5(A))

static double scal_prod(const double *x, const double *y, const int size)
{
int nblock = size/10;
int remainder = size - nblock*10;
double sum = 0;
for (int i=0; i<nblock; ++i) 
  {
  DO10( sum+= (*x)*(*y); ++x; ++y;)
  }
 for (int i=0; i<remainder; ++i) {sum+= (*x)*(*y); ++x; ++y;}
return sum;
}  


static double three_scal_prod(double *x, double *y, double *z, const int size)
{
int nblock = size/10;
int remainder = size - nblock*10;
double sum = 0;
for (int i=0; i<nblock; ++i) 
  {
    DO10( sum+= (*x)*(*y)*(*z); ++x; ++y; ++z)
  }
 for (int i=0; i<remainder; ++i) {sum+= (*x)*(*y)*(*z); ++x; ++y;++z;}
return sum;
}  


int KernelFit::FitDifferentialBackground(ImagePair &ImPair, const double NSig)
{  
   
  if (optParams.SepBackVar.Degree == -1) return 0;

  int nterms = optParams.SepBackVar.Nterms();

  diffbackground.resize(nterms);

  Mat A(nterms,nterms);
  Vect B(nterms);
  Vect monom(nterms);

  Pixel bestMean,bestSig,worstMean,worstSig;
  const Image &bestImage = ImPair.BestImage();
  const Image &worstImage = ImPair.WorstImage();
  const Frame &dataFrame = ImPair.CommonFrame();
  bestImage.SkyLevel(dataFrame,&bestMean, &bestSig);
  worstImage.SkyLevel(dataFrame,&worstMean, &worstSig);
  
  double cut1  = NSig*bestSig;
  double cut2  = NSig*worstSig;

  int ibeg,iend,jbeg,jend;
  ibeg = int(dataFrame.xMin);
  iend = int(dataFrame.xMax);
  jbeg = int(dataFrame.yMin);
  jend = int(dataFrame.yMax);

  int jump = 10;
  for (int j=jbeg ; j< jend ; j+= jump)
  for (int i=ibeg ; i< iend ; i+= jump)
    {
      Pixel p1 = worstImage(i,j);
      if (fabs(p1-bestMean) > cut1) continue;
      Pixel p2 = bestImage(i,j);
      if (fabs(p2-worstMean) > cut2) continue;
      for (int q1=0; q1<nterms; ++q1)
	{
	  monom(q1) = optParams.SepBackVar.Value(double(i), double(j),q1);

	  for (int q2 = q1; q2<nterms; ++q2) A(q1,q2) += monom(q1)*monom(q2);
	  B(q1) += monom(q1)*(p1-p2);
	}
    }
  /* symetrize */
  for (int q1=0; q1<nterms; ++q1) 
    for (int q2 = q1+1; q2<nterms; ++q2) A(q2,q1) = A(q1,q2); 
  if (cholesky_solve(A,B,"U")!= 0)
    {
      cerr << " could not compute differential background !!!!" << endl;
      return 0;
    }
  
  cout << setprecision(10);
  cout << " separately fitted differential background " << endl;
  cout << " ----------------------------------------- " << endl;
  for (int q1=0; q1< nterms; ++q1) cout << B(q1) << " " ;
  cout << endl;

  for (int q1=0; q1<nterms; ++q1)
    diffbackground[q1] = B(q1);
  return 1;
}



#define OPTIMIZED /* means pushing pointers by hand ... */

// #define DO3(I) I;I;I;
// #define DO9(I) DO3(DO3(I))
// #define DO19(I) DO9(I);DO9(I);I

#include <time.h>


//! convolves the best image with the current kernel (Usually set by DoTheFit).
/*! The UpdateKern parameter is the pixel range over which the kernel will not
be updated.*/

void KernelFit::ImageConvolve(const Image &Source, Image &Out, int UpdateKernStep)
{
  if (!fitDone)
    throw(PolokaException("KernelFit::ImageConvolve :  the kernel fit was not done yet ")); 
  clock_t tstart = clock();
  assert(&Source != &Out);

  if (solution.size()==0)
    {
      cerr << "BestImageConvolve : the kernel fit was not done yet... or was impossible .. no convolution" << endl;
      return;
    }
  // copy the input image so that side bands are filled with input values.
  Out = Source; 

  int ksx = optParams.HKernelSize;
  int ksy = optParams.HKernelSize;
  int startx = ksx;
  int starty = ksy;
  int endx = Source.Nx() - ksx;
  int endy = Source.Ny() - ksy;
  Kernel kern(ksx,ksy);
  int npix = kern.Nx();
  int nxregions = (UpdateKernStep) ? (Source.Nx()/UpdateKernStep)+1 : 1;
  int nyregions = (UpdateKernStep) ? (Source.Ny()/UpdateKernStep)+1 : 1;
  cout << " Convolving best image" << endl;
  // some printout:
  KernCompute(kern, Source.Nx()/2, Source.Ny()/2);
  {cout <<  " kernel caracteristics at i j " <<  Source.Nx()/2 << ' ' << Source.Ny()/2; kern.dump_info();}
  for (int i = 1 ; i <= 3 ; i++)
    for (int j = 1 ; j <= 3 ; j++)
      {
	KernCompute(kern, i*Source.Nx()/4, j*Source.Ny()/4);
	cout <<  " kernel caracteristics at i j " <<  i*Source.Nx()/4 
	     << ' ' << j*Source.Ny()/4; kern.dump_info();
      }
  
  // on the road again
  for (int iry = 0; iry < nyregions; ++iry)
    for (int irx = 0; irx < nxregions; ++irx)
      {
	int sj = max(iry*UpdateKernStep,starty);
	int ej = min(sj+UpdateKernStep,endy);
	int si = max(irx*UpdateKernStep, startx);
	int ei = min(si+UpdateKernStep, endx);
	double xc = (si+ei)/2;
	double yc = (sj+ej)/2;
	KernCompute(kern, xc, yc);
	//	double kernSum = kern.sum();
	
	if (((irx == 0) || (irx == nxregions-1)) &&  ((iry ==0) || (iry == nyregions - 1)))
	  {cout <<   " kernel caracteristics at i j " << xc << ' ' << yc; kern.dump_info();}
	for (int j = sj; j < ej; ++j)
	  for (int i = si; i < ei; ++i)
	    {
	      double sum = 0;
#ifndef OPTIMIZED
	      for (int jk = -ksy; jk <=ksy; ++jk)
		for (int ik = -ksx; ik <=ksx; ++ik)
		  sum += kern(ik,jk)*Source(i-ik,j-jk);
#else
	      DPixel *pk = kern.begin();
	      for (int jk = -ksy; jk <= ksy; ++jk)
		{
		  Pixel *ps = &Source(i+ksx, j-jk);
		  for (int k=npix; k; --k) {sum += (*pk) * (*ps); ++pk ; --ps;}
		}
#endif
	      // A REGARDER:
	      //      sum -= kernSum*BestImageBack;
  //  if (kern.begin() + kern.Nx()*kern.Ny() - pk ) cout << " BestImageImageConvol catastrophe...." << endl;
      Out(i,j) = sum;
    }
  }
/* account for differential background. We do it for the whole image, including
   side bands where the convolution did not go but where anyway initialized
   to the input value. TODO : correct side bands for photometric ratio */

 int sx = Out.Nx();
 int sy = Out.Ny();

 clock_t tend = clock();
 cout << " CPU for convolution " << float(tend-tstart)/float(CLOCKS_PER_SEC) << endl;


 if (optParams.BackVar.Nterms()==0) // sep fit background
   for (int j=0; j < sy; ++j) for (int i=0; i < sx ; ++i) Out(i,j) += SepBackValue(i,j);
 else
   for (int j=0; j < sy; ++j) for (int i=0; i < sx ; ++i) Out(i,j) += BackValue(i,j);

 tend = clock();
 cout << " CPU for convolution+back sub " << float(tend-tstart)/float(CLOCKS_PER_SEC) << endl;
}

void KernelFit::VarianceConvolve(const Image &Source, Image&Out, int UpdateKernStep)
{
  if (!fitDone)
    throw(PolokaException("KernelFit::ImageConvolve :  the kernel fit was not done yet ")); 
  assert(&Source != &Out);
  clock_t tstart = clock();
if (solution.size()==0)
   {
     cerr << "Variance Convolve : the kernel fit was not done yet... or was impossible .. no convolution" << endl;
     return;
   }
int ksx = optParams.HKernelSize;
int ksy = optParams.HKernelSize;
int startx = ksx;
int starty = ksy;
int endx = Source.Nx() - ksx;
int endy = Source.Ny() - ksy;
Kernel kern(ksx,ksy);
int npix = kern.Nx();
int nxregions = (UpdateKernStep) ? (Source.Nx()/UpdateKernStep)+1 : 1;
int nyregions = (UpdateKernStep) ? (Source.Ny()/UpdateKernStep)+1 : 1;
cout << " Convolving variance" << endl;
// some printout:
KernCompute(kern, Source.Nx()/2, Source.Ny()/2);
 {cout <<  " kernel caracteristics at i j " <<  Source.Nx()/2 << ' ' << Source.Ny()/2; kern.dump_info();}
 for (int i = 1 ; i <= 3 ; i++)
   for (int j = 1 ; j <= 3 ; j++)
     {
       KernCompute(kern, i*Source.Nx()/4, j*Source.Ny()/4);
       kern *= kern;
       cout <<  " kernel^2 caracteristics at i j " <<  i*Source.Nx()/4 
	    << ' ' << j*Source.Ny()/4; kern.dump_info();
     }

 // on the road again
for (int iry = 0; iry < nyregions; ++iry)
for (int irx = 0; irx < nxregions; ++irx)
  {
  int sj = max(iry*UpdateKernStep,starty);
  int ej = min(sj+UpdateKernStep,endy);
  int si = max(irx*UpdateKernStep, startx);
  int ei = min(si+UpdateKernStep, endx);
  double xc = (si+ei)/2;
  double yc = (sj+ej)/2;
  KernCompute(kern, xc, yc);
  
  // This is a variance :
  kern *= kern;

  for (int j = sj; j < ej; ++j)
  for (int i = si; i < ei; ++i)
    {
      double sum = 0;
#ifndef OPTIMIZED
      for (int jk = -ksy; jk <=ksy; ++jk)
      for (int ik = -ksx; ik <=ksx; ++ik)
         sum += kern(ik,jk)*Source(i-ik,j-jk);
#else
      DPixel *pk = kern.begin();
      for (int jk = -ksy; jk <= ksy; ++jk)
	{
	  Pixel *ps = &Source(i+ksx, j-jk);
 	  for (int k=npix; k; --k) {sum += (*pk) * (*ps); ++pk ; --ps;}
	}
#endif
      Out(i,j) = sum;
    }
  }
clock_t tend = clock();
cout << " CPU for convolution " << float(tend-tstart)/float(CLOCKS_PER_SEC) << endl;
}


static double my_pow(const double X, const int k)
{
  if (k==0) return 1;
  double res = X;
  for (int i=1; i<k; ++i) res *= X;
  return res;
}


static void PolGaussKern(Kernel &Kern, double Sig, int DegX, int DegY)
{
  /* fills the kernel array with values : exp (-(x**2+y**2)/(2.sigma**2))*(x**degx)*(y**degy) */
int s = Kern.HSizeX();
int totSize = 2*s+1;
double *vx = new double[totSize];
double *vy = new double[totSize];
double alpha = 1./(2.*Sig*Sig);
for (int i=0; i<totSize; ++i)
  {
  double x = double(i-s);
  double valg = exp(-x*x*alpha);
  vx[i] = valg*my_pow(x,DegX);
  vy[i] = valg*my_pow(x,DegY);
  }

DPixel *p = Kern.begin();
for (int j=0;j<totSize;++j)
  for (int i=0; i<totSize; ++i)
    {*p = vx[i]*vy[j]; ++p;}
//#warning on ne peut pas relire des kernels avec ce code.
// Kern *= 1/(Kern.sum2());

delete [] vx;
delete [] vy;
}
 

static void SetDelta(Kernel &Kern)
{
Kern.Zero();
Kern(0,0) = 1.0;
}


static vector<int> toto;




void KernelFit::BasisFill(Mat *BasisTransfo)
{
  if (Basis.size() != 0)
    {
      cout << "entering BasisFill with a non empty Basis array ... " << endl;
      Basis.clear(); 
    }
  /* all kernels should have the same size */
  Basis.push_back(Kernel(optParams.HKernelSize, optParams.HKernelSize));
  SetDelta(Basis[0]);
  int oldprec = cout.precision();
  cout << setprecision(10);
  cout << "  kernels parameters (sigma, degree) : ";
  for (int iwidth = 0; iwidth < optParams.NGauss ; ++iwidth) 
    cout << '(' <<optParams.Sigmas[iwidth] << ',' << optParams.Degrees[iwidth] << ')';
  cout << endl;  
  
  for (int iwidth = 0; iwidth < optParams.NGauss ; ++iwidth)
    {
      int maxdeg = optParams.Degrees[iwidth];
      /* do not change the ordering here without changing the IO 
	 version number and code the handling of "old" saved kernels */
      for (int degx = 0; degx <= maxdeg; ++degx)
	for (int degy =0; degy<= maxdeg-degx; ++degy)
	  {
	    Kernel kern(optParams.HKernelSize, optParams.HKernelSize);
	    //        cout << " degx degy " << degx << ' ' << degy << endl;
	    PolGaussKern(kern,optParams.Sigmas[iwidth], degx, degy);
	    Basis.push_back(kern);
	  }
    }
  mSize = Basis.size() * optParams.KernVar.Nterms() + optParams.BackVar.Nterms();
  cout << " BasisFill : Number of basic kernels : " << Basis.size()
       << " number of fitted coefficients "  << mSize <<endl;
  
  
  if (optParams.OrthogonalBasis)
    {
      unsigned n = Basis.size();
      int kernelLength=Basis[0].Nx() * Basis[0].Ny();
      vector<Kernel> oldBasis(Basis);
      Basis.clear();
      Basis.reserve(n);
      for (unsigned i=0; i<n; ++i)
	{
	  Kernel kern(oldBasis[i]);
	  for (unsigned j=0; j<i ; ++j)
	    {
	      double coeff = scal_prod(kern.begin(), Basis.at(j).begin(), kernelLength);
	      kern -= Basis.at(j)*coeff;
	    }	      
	  double norm = kern.sum2();
	  assert(norm>0);
	  kern *= 1./sqrt(norm);
	  Basis.push_back(kern);
	}
      if (BasisTransfo)
	{
	  BasisTransfo->allocate(Basis.size(), Basis.size());
	  for (unsigned j=0;j<Basis.size(); ++j)
	    for (unsigned i=0; i<Basis.size() ; ++i)
	      (*BasisTransfo)(i,j) = scal_prod(oldBasis[i].begin(), Basis[j].begin(), kernelLength);
	}	  
    } // end if (optParams.OrthogonalBasis)
  else
    if (BasisTransfo) // allocate and fill identity matrix
      {
	BasisTransfo->allocate(Basis.size(), Basis.size());
	for (unsigned j=0;j<Basis.size(); ++j) (*BasisTransfo)(j,j) = 1;
      }
  cout << setprecision(oldprec);
}

/* this routine is only useful to study the fit */
void KernelFit::ParameterGroups(Mat &Groups) const
{
  int ngroups = 1; // delta
  int spatVar = optParams.KernVar.Nterms();
  // element basis
  for (int iwidth = 0; iwidth < optParams.NGauss ; ++iwidth)
    {
      int maxdeg = optParams.Degrees[iwidth];
      ngroups += maxdeg+1;
    }
  // spatial variation
  ngroups += optParams.KernVar.Degree;

  Groups.allocate(ngroups,mSize);
  // now fill  
  // delta
  for (int is = 0; is< spatVar; ++is) Groups(0,KernIndex(0,is) ) = 1;
  int igroup = 1;
  // Gaussians*pol
  int start = 1; // skip delta
  for (int iwidth = 0; iwidth < optParams.NGauss ; ++iwidth)
    {
      int maxdeg = optParams.Degrees[iwidth];
      int end = start+((maxdeg+1)*(maxdeg+2))/2;
      int is = start;
      for (int deg=0; deg <= maxdeg; deg++)
	{
	  for (int ik = is; ik<end; ++ik)
	    for (int is = 0; is< spatVar; ++is)
	      Groups(igroup, KernIndex(ik,is) ) = 1;
	  is += (deg+1);
	  igroup++;
	}
      is = start;
    }
  // spatial variation
  for (int isDeg=1; isDeg<= optParams.KernVar.Degree; isDeg++)
    {
      for (unsigned ik=0; ik < Basis.size(); ++ik)
	for (int is=0; is<spatVar; ++is)
	  {
	    int xdeg = optParams.KernVar.Xdeg[is];
	    int ydeg = optParams.KernVar.Ydeg[is];
	    if (xdeg+ydeg>= isDeg) Groups(igroup, KernIndex(ik,is)) = 1;
	  }
      igroup++;
    }
}
  

vector<double> DeltaChi2(const Mat &A, const Vect &B, const Mat &Groups)
{
  Vect Sol(B);
  Mat Am1(A);
  cholesky_solve(Am1,Sol, "U");
  cholesky_invert(Am1,"U");
  Am1.Symmetrize("U");
  vector<double> res;
  for (unsigned ig=0; ig<Groups.SizeX(); ++ig)
    {
      vector<unsigned> indices;
      for (unsigned j=0; j<Groups.SizeY(); ++j) if (Groups(ig,j) != 0) indices.push_back(j);
      unsigned dSize = indices.size();
      Mat D(dSize,dSize);
      Vect E(dSize);
      for (unsigned i=0; i<dSize; ++i)
	{
	  for (unsigned j=0; j<dSize; ++j)
	    D(i,j) = Am1(indices[i], indices[j]);
	  E(i) = Sol(indices[i]);
	}
      Vect Ep(E);
      if (cholesky_solve(D,E,"U") != 0)
	{
	  cout << "aie aie aie " << endl;
	  continue;
	}
      res.push_back(Ep*E);
    }
  return res;
}

static double image_scal_prod(const DImage &Vignette, const DImage &w, const Image& I, int xs, int ys)
{
double sum = 0;
int xsize = Vignette.Nx();
int ysize = Vignette.Ny();
for (int j=0; j< ysize; ++j)
  {
  for (int i=0; i< xsize; ++i)
    {     
      sum += double(I(i+xs,j+ys)) * Vignette(i,j)*w(i,j);
    } 
  }
return sum;
}

#ifdef STORAGE
static double image_sum(const Image* I, const double IBack, int xs, int ys, int xe, int ye)
{
double sum = 0;
for (int j=ys; j<=ye ; ++j)
  {
  for (int i=xs; i<=xe; ++i)
    {
    sum += double((*I)(i,j)) - IBack;
    } 
  }
return sum;
}
#endif

class FitWorkSpace {

/* allocate vectors of DImages to store various things needed for the fit 
   (e.g. a_stamp x Basis[i]) */
 public :
  vector<DImage> convolutions;
  vector<DImage> backStamps;

  FitWorkSpace(const KernelFit&);
  
  // no need for an explicit destructor : since ~vector will do the job.
  
};


FitWorkSpace::FitWorkSpace(const KernelFit& Fit)
{
  int nkern = Fit.Basis.size();
  convolutions.assign(nkern,DImage());
  int convolvedSize = Fit.optParams.ConvolvedSize();
  for (int i=0; i<nkern; ++i) 
    convolutions[i].Allocate(convolvedSize,convolvedSize);
  int nBackStamps = Fit.optParams.BackVar.Nterms();
  backStamps.assign(nBackStamps, DImage());
  for (int i=0; i<nBackStamps; ++i) 
    backStamps[i].Allocate(convolvedSize,convolvedSize);
}

static FitWorkSpace *work = NULL;



void KernelFit::OneStampNoiseMatrix(const Stamp &AStamp, 
				    const Image &WorstImage, Mat &M) const
{
  
  int nkern = Basis.size();
  DImage  bv_conv_w;
  unsigned cs = AStamp.bestVariance.Nx()-AStamp.weight.Nx()+1;
  bv_conv_w.Allocate(cs,cs);
  DImage mirroredVariance;
  Mirror(AStamp.bestVariance, mirroredVariance);
  Convolve(bv_conv_w, mirroredVariance, AStamp.weight);
  // DEBUG
  if (!FileExists("bestVar.fits"))
    {
      AStamp.bestVariance.writeFits("bestVar.fits");
      AStamp.weight.writeFits("weight.fits");
      bv_conv_w.writeFits("bv_conv_w.fits");
    }

  /* values of monomials such as (xc**n)*(yc**m), for this stamp  : */
  // for the smooth kernel variation
  double *spatialCoeff = new double [optParams.KernVar.Nterms()];
  for (unsigned int i = 0; i  < optParams.KernVar.Nterms(); ++i) 
    spatialCoeff[i] = optParams.KernVar.Value(AStamp.xc, AStamp.yc, i);

  /* compute here sum_pix ((Ki^2 \conv Var(Best))* weight),
     which turns out to be equal to : sum_pix ((Ki^2)* (Flip(Var(Best)) \conv weight)) */
  int nPixKern = (nkern >0) ? Basis[0].Nx()*Basis[0].Ny() : 0;
  for (int ik=0; ik < nkern; ++ik)
    {
      DImage product(Basis[ik]);
      product *=bv_conv_w;
      for (int jk=0; jk <=ik; ++jk)
	{
	  double integral = scal_prod(Basis[jk].begin(), product.begin(), nPixKern);
	  for (unsigned int is =0; is < optParams.KernVar.Nterms(); ++is)
	    {
	      int im = KernIndex(ik,is);
	      for (unsigned int js =0; js < optParams.KernVar.Nterms(); ++js)
		{
		  int jm = KernIndex(jk,js);
		  double val = integral*spatialCoeff[is]*spatialCoeff[js];
		  M(im,jm) += val; 
		}
	    }
	}
    }
  delete [] spatialCoeff;
  //symmetrise
  for (size_t i=0; i<mSize; ++i) 
  for (size_t j=i+1; j<mSize; ++j) 
    M(i,j) = M(j,i);


}

void KernelFit::OneStampMAndB(const Stamp &AStamp, const Image &WorstImage, Mat &M, Vect &B)
{
  /* computes the contributions of AStamp to the linear system
     and increments accordingly M and B, which are *NOT* zeroed here.
     Note that only the "U" part of M is used on input, 
     and M is symetrized on output, so that the UorL choice remains confined
     to this routine.
  */
  // useless : no printouts in the routine .
  int oldprec = cout.precision();
 cout << setprecision(10);
 ios::fmtflags  old_flags = cout.flags(); 
 cout << resetiosflags(ios::fixed) ;
 cout << setiosflags(ios::scientific) ;

 int nkern = Basis.size();
/* assume that all kernels have the same size: */
 //int ksize = optParams.HKernelSize;
/* assume that all stamps have the same size: */
 //int stampsize = AStamp.Nx();
 int convolvedSize = optParams.ConvolvedSize(); 

 /* 'convolutions' and backStamps are parts of "work" and should be allocated somewhere else
    so that this is allocated once for all stamps :*/

 int hConvolvedSize = convolvedSize/2;
 int convolvedPix = convolvedSize*convolvedSize;
 //DEBUG
 assert(AStamp.weight.Nx() == convolvedSize);

 double *spatialCoeff = new double [optParams.KernVar.Nterms()];
 // double *backCoeff = new double [optParams.BackVar.Nterms()];


 /* values of monomials such as (xc**n)*(yc**m), for this stamp  : */
 // for the smooth kernel variation
 for (unsigned int i = 0; i  < optParams.KernVar.Nterms(); ++i) 
     spatialCoeff[i] = optParams.KernVar.Value(AStamp.xc, AStamp.yc, i);
 //for the background, actually compute monomials for all stamp pixels:
 for (unsigned int ib =0; ib < optParams.BackVar.Nterms(); ++ib) 
   {
     DImage &backStamp = work->backStamps[ib];
     for (int j= 0; j < convolvedSize ; ++j)
     for (int i= 0; i < convolvedSize ; ++i)
       {
	 int xi = i+ AStamp.xc - hConvolvedSize;
	 int yj = j+ AStamp.yc - hConvolvedSize;
	 backStamp(i, j) = optParams.BackVar.Value(xi, yj, ib);
       }
   }

    /* contributions to the matrix m */
    /* background-background terms */
 for (unsigned int ib = 0; ib < optParams.BackVar.Nterms(); ++ib)
   {
   for (unsigned int jb =0; jb<=ib; ++jb)   
     {
       /* BackIndex is an increasing function of its argument, 
	  so we are filling the "U" part (j<=i) */
       M(BackIndex(ib),BackIndex(jb)) += 
       three_scal_prod(work->backStamps[ib].begin(), work->backStamps[jb].begin(), 
		       AStamp.weight.begin(), convolvedPix);
     }
   }


    /* kernel - kernel terms */
 for (int ik=0; ik < nkern; ++ik)
   {
   Convolve(work->convolutions[ik], AStamp.bestPixels, Basis[ik]);
   for (int jk=0; jk <=ik; ++jk)
     {
     double integral = three_scal_prod(work->convolutions[ik].begin(), 
				       work->convolutions[jk].begin(),
				       AStamp.weight.begin(),
				       convolvedPix);
     for (unsigned int is =0; is < optParams.KernVar.Nterms(); ++is)
       {
       int im = KernIndex(ik,is);
       for (unsigned int js =0; js < optParams.KernVar.Nterms(); ++js)
	 {
	 int jm = KernIndex(jk,js);
         M(im,jm) += integral*spatialCoeff[is]*spatialCoeff[js];  // alard eq (3)
         }
       }
     }
   }

 for (int ik=0; ik < nkern; ++ik)
      /* kernel-background terms */
   for (unsigned int is = 0; is < optParams.KernVar.Nterms(); ++is)
     {
     int im = KernIndex(ik,is);
     for (unsigned int jb=0; jb < optParams.BackVar.Nterms(); ++jb) 
       {
       int jm = BackIndex(jb);
       /* fill only M(i,j) for j<=i,
	  and background is placed after kernel in params */
       M(jm,im) += spatialCoeff[is]*
	 three_scal_prod(work->convolutions[ik].begin(), 
			 work->backStamps[jb].begin(), 
			 AStamp.weight.begin(),  
			 convolvedPix);
       }
     }

      /* compute contributions to b */
 for (int ik=0; ik < nkern; ++ik)
   {
      /* kernel terms */
   double bintegral = image_scal_prod(work->convolutions[ik], AStamp.weight, WorstImage, 
				      AStamp.xc - hConvolvedSize, AStamp.yc - hConvolvedSize);
   for (unsigned int is =0; is < optParams.KernVar.Nterms(); ++is)
     {
     int im = KernIndex(ik,is);
     B(im) += bintegral * spatialCoeff[is]; // alard eq (4)
     }
   } /* end of for (ik */

      /* background terms of b */
 for (unsigned int ib=0; ib < optParams.BackVar.Nterms(); ++ib)
   {
     B(BackIndex(ib)) += 
     image_scal_prod(work->backStamps[ib], AStamp.weight, WorstImage,
		     AStamp.xc - hConvolvedSize, AStamp.yc - hConvolvedSize);
   }
 delete [] spatialCoeff;
 // delete [] backCoeff;
 /* return a symetrized matrix ... safer than assuming anything 
    about the calling routine*/
for (size_t i=0; i<mSize; ++i) 
  for (size_t j=i+1; j<mSize; ++j) 
    M(i,j) = M(j,i);

 cout << setprecision(oldprec);
 cout.flags(old_flags);
}





static void KernLinComb(Kernel &Result, const vector<Kernel> &VK, double *Coeffs)
{
int size = Result.Nx() * Result.Ny();
int nkern = int(VK.size());
Result.Zero();
for (int ik=0; ik < nkern; ++ik)
  {
  DPixel *k = VK[ik].begin();
  DPixel *r = Result.begin();
  double coeff = Coeffs[ik];
  for (int i=0; i< size; ++i) {*r += *k * coeff; ++r ; ++k;}
  }
}


#include "toadscards.h"
#include "datacards.h"
#include "fileutils.h" // for FileExists

OptParams::OptParams()
{
NGauss = 3;
Sigmas.resize(NGauss);
Degrees.resize(NGauss);
 Sigmas[0] = 0.7; Sigmas[1] = 1.5; Sigmas[2] = 2.;
Degrees[0] = 6; Degrees[1] = 4; Degrees[2] = 2;
//Degrees[0] = 1; Degrees[1] = 4; Degrees[2] = 2;
HStampSize = 15;
HKernelSize = 9;
MaxStamps = 1000;
KernVar.SetDegree(2); // degree of spatial variations of the kernel 
BackVar.SetDegree(2); //degree of spatial variations of the background
SepBackVar.SetDegree(-1); //degree of spatial variations of the background if you want to fit it separately
 NSig = 4;
UniformPhotomRatio = true;
OrthogonalBasis = true;
 SubtractNoise = true;
 string dataCardsName = DefaultDatacards();
 
 if (FileExists(dataCardsName))
   {
     cout << " KernelFit uses datacards : " << dataCardsName << endl;
     DataCards cards(dataCardsName);
     if (cards.HasKey("KFIT_MAX_STAMPS")) {
       MaxStamps = cards.IParam("KFIT_MAX_STAMPS");
     }
     if (cards.HasKey("KFIT_SIG_GAUSS"))
       {
	 int nGauss = cards.NbParam("KFIT_SIG_GAUSS");
	 int nGaussDeg = cards.NbParam("KFIT_DEG_GAUSS");
	 if (nGaussDeg != nGauss)
	   {
	     cerr << dataCardsName 
		  << " : don't find the same number of items in KFIT_SIG_GAUSS and KFIT_DEG_GAUSS" 
		  << endl;
	     cerr << " cutting the longest one  to the size of the shortest one "<< endl;
	     nGauss = min(nGauss,nGaussDeg);
	   }
	 if (nGauss > 0) 
	   {
	     NGauss = nGauss;
	     Sigmas.resize(NGauss);
	     Degrees.resize(NGauss);
	     for (int i=0; i<NGauss; ++i)
	       {
		 Sigmas[i] = cards.DParam("KFIT_SIG_GAUSS", i);
		 Degrees[i] = cards.IParam("KFIT_DEG_GAUSS",i);
	       }	   
	   }
       } // HasKey("KFIT_SIG_GAUSS")
     if (cards.HasKey("KFIT_NSIG")) NSig = cards.DParam("KFIT_NSIG");
     if (cards.HasKey("KFIT_KERNVAR_DEG")) 
       {
	 int deg = cards.IParam("KFIT_KERNVAR_DEG");
	 KernVar.SetDegree(deg); // degree of spatial variations of the kernel 
       }
     if (cards.HasKey("KFIT_BACKVAR_DEG")) 
       {
	 int deg = cards.IParam("KFIT_BACKVAR_DEG");
	 BackVar.SetDegree(deg); // degree of spatial variations of the background
       }
     if (cards.HasKey("KFIT_SEPBACKVAR_DEG")) 
       {
	 int deg = cards.IParam("KFIT_SEPBACKVAR_DEG");
	 SepBackVar.SetDegree(deg); // degree of spatial variations of the background
       }
     if (cards.HasKey("KFIT_UNIFORM_PHOTOM_RATIO")) 
       UniformPhotomRatio = cards.IParam("KFIT_UNIFORM_PHOTOM_RATIO") == 1;
     if (cards.HasKey("KFIT_SUBTRACT_NOISE")) 
       SubtractNoise = cards.IParam("KFIT_SUBTRACT_NOISE") == 1;

   } // if (has datacards)
}
  

void OptParams::OptimizeSizes(double BestSeeing, double WorstSeeing)
{
 int oldprec = cout.precision();
 cout << setprecision(10);
cout << " Choosing sizes with bestseeing = " <<  BestSeeing << " WorstSeeing = " << WorstSeeing << endl;
if ( BestSeeing > WorstSeeing) swap(BestSeeing, WorstSeeing);
double kernSig = max(sqrt(WorstSeeing*WorstSeeing - BestSeeing*BestSeeing),0.4);
 cout << " Expected kernel sigma  : " << kernSig << endl ; 
HKernelSize = max(int( ceil (NSig * kernSig)),4);
 for (int i=0; i<NGauss; ++i) Sigmas[i] *= kernSig;
/* convolvedStampSize is the size of the stamps that enter the chi2. 2 is added
to get some area to estimate the differential background. */
int convolvedStampSize = int(ceil(NSig * WorstSeeing)) + 2;
HStampSize = HKernelSize + convolvedStampSize;
 cout << setprecision(oldprec);
}

void OptParams::OptimizeSpatialVariations(const int NumberOfStars)
{
  int degree = KernVar.Degree;
  if (NumberOfStars < 75) degree--;
  if (NumberOfStars < 50) degree--;
  // only lower greees, not raise it from datacards value
  if (KernVar.Degree > degree) 
    {
      cout << " Lowering spatial kernel (deg = " << degree 
	   << ") variations with NumberOfStars = " << NumberOfStars << endl;
      KernVar.SetDegree(degree);
    }
  //  if (BackVar.Degree > degree) BackVar.SetDegree(degree); not so useful

  if (NumberOfStars < 25)
    {
      cout << " Choosing spatial variations with NumberOfStars = " 
	   << NumberOfStars << endl;
      cout << " lowered degree of polynomials associated to gaussians to ";
      for (int i=0; i<NGauss; ++i) {Degrees[i] /= 2; cout << Degrees[i] << ' ';}
      cout << endl;
    }

  if (NumberOfStars < 8) 
    {
      cout << " lowered degree of polynomials associated to gaussians to ";
      for (int i=0; i<NGauss; ++i) {Degrees[i] /= 2; cout << Degrees[i] << ' ';}
      cout << endl;
    }
}


void OptParams::read(istream& stream)
{
  string tmp_str;
  int version;
  stream >> tmp_str >> version;
  if (version>1) throw(PolokaException(" unknown version number in OptParams::read"));
  read_member(stream, tmp_str, HKernelSize);
  read_member(stream, tmp_str, NGauss);
  read_member(stream, tmp_str, Sigmas);
  read_member(stream, tmp_str, Degrees);
  read_member(stream, tmp_str, NSig);
  read_member(stream, tmp_str, KernVar);
  read_member(stream, tmp_str, BackVar);
  read_member(stream, tmp_str, SepBackVar);
  read_member(stream, tmp_str, HStampSize);
  read_member(stream, tmp_str, MaxStamps);
  read_member(stream, tmp_str, UniformPhotomRatio);
  if (version >=1) read_member(stream, tmp_str, OrthogonalBasis);
  else OrthogonalBasis=false;
}


void OptParams::write(ostream& stream) const
{
  static int OptParams_version = 1;
  stream << "[OptParams] " <<   OptParams_version ;
  write_member(stream, "HKernelSize", HKernelSize);
  write_member(stream, "NGauss", NGauss);
  write_member(stream, "Sigmas", Sigmas);
  write_member(stream, "Degrees", Degrees);
  write_member(stream, "NSig", NSig);
  write_member(stream, "KernVar", KernVar);
  write_member(stream, "BackVar", BackVar);
  write_member(stream, "SepBackVar", SepBackVar);
  write_member(stream, "HStampSize", HStampSize);
  write_member(stream, "MaxStamps", MaxStamps);
  write_member(stream, "UniformPhotomRatio", UniformPhotomRatio);
  write_member(stream, "OrthogonalBasis", OrthogonalBasis);
}




void XYPower::SetDegree(const int DegreeValue)
{
Degree = DegreeValue;
nterms = (Degree+1)*(Degree+2)/2;
Xdeg.resize(nterms);
Ydeg.resize(nterms);
int q = 0;
for (int xdeg=0; xdeg <=Degree; ++xdeg)
  for (int ydeg = 0; ydeg <= Degree - xdeg; ++ydeg)
    {
    Xdeg[q] = xdeg;
    Ydeg[q] = ydeg;
    ++q;
    }
}


double XYPower::Value(const double X, const double Y, const unsigned q) const
{
  if ((unsigned int)q>=Nterms()) {cerr << "  XYPower::Value ..."  << endl; abort();}
// my_pow(double,int) is about 5 times faster than pow(double,doublke)
  return my_pow(X/100.,Xdeg[q])*my_pow(Y/100.,Ydeg[q]); 
}


void XYPower::ApplyToImage(Image&I, double Factor, const vector<double> &ParamVal) const
{
  assert(ParamVal.size() >= nterms);
  double *xpow = new double[I.Nx()*nterms];
  double *ypow = new double[I.Ny()*nterms];
  for (int i=0; i<I.Nx(); ++i)
    for (int k=0; k<nterms; ++k) xpow[i*nterms+k] = my_pow(double(i)/100.,Xdeg[k])*Factor*ParamVal[k];
  for (int j=0; j<I.Ny(); ++j)
    for (int k=0; k<nterms; ++k) ypow[j*nterms+k] = my_pow(double(j)/100.,Ydeg[k]);
  for (unsigned j=0; j<I.Ny(); ++j)
    for (unsigned i=0; i<I.Nx(); ++i) 
      {
	double val = 0;
	double *px = &xpow[i*nterms];
	double *py = &ypow[j*nterms];
	for (unsigned k=0; k<nterms; ++k) {val += (*px)*(*py); ++px; ++py;}
	  I(i,j) += val;
      }
  delete [] xpow; delete [] ypow;
}





void XYPower::read(istream& stream)
{
  string tmp_str;
  int version;
  stream >> tmp_str >> version;
  read_member(stream, tmp_str, Degree);
  read_member(stream, tmp_str, Xdeg);
  read_member(stream, tmp_str, Ydeg);
  nterms = Xdeg.size();
}

void XYPower::write(ostream& stream) const
{
  stream << "[XYPower] " << 0;
  write_member(stream, "Degree", Degree);
  write_member(stream, "Xdeg", Xdeg);
  write_member(stream, "Ydeg", Ydeg);
}





void KernelFit::KernCompute(Kernel &Kern, const double X, const double Y) const
{
unsigned int nKern = Basis.size();
double *coeff = new double[nKern];
for (unsigned int i=0; i<nKern; ++i) coeff[i] = 0;
for (unsigned int is =0; is < optParams.KernVar.Nterms(); ++is)
  {
  double spatialCoeff = optParams.KernVar.Value(X,Y,is);
  for (unsigned int ik=0; ik< nKern; ++ik)
    {  
      coeff[ik] +=  solution.at(KernIndex(ik, is))*spatialCoeff;
    }
  }
KernLinComb(Kern, Basis, coeff);
delete [] coeff;
}

void KernelFit::KernAllocateAndCompute(Kernel &Kern, const double X, const double Y) const
{
  Kern.allocate(HKernelSizeX(), HKernelSizeY());
  KernCompute(Kern, X, Y);
}

static void mean_median_sigma(double *values, const int nval, double &mean,  double &median, double &sigma)
{
mean =0;
sigma =0;
median = DArrayMedian(values, nval);
for (int i=0; i<nval-1 ; ++i)
  {
  mean += values[i];
  sigma += values[i]*values[i];
  }
mean /= double(nval);
sigma = sigma/double(nval) - mean*mean;
if (sigma>0)  sigma = sqrt(sigma); else sigma = 0;
}



double KernelFit::StampChi2(Stamp &stamp, const Image &WorstImage)
{
Kernel kern(optParams.HKernelSize, optParams.HKernelSize);
int  xc = stamp.xc;
int  yc = stamp.yc;
KernCompute(kern, xc, yc);
int convolvedSize = optParams.ConvolvedSize(); /* should be odd */
int hConvolvedSize = convolvedSize/2;
DImage R_conv_K(convolvedSize,convolvedSize);
Convolve(R_conv_K, stamp.bestPixels, kern);
double chi2 = 0;
int xs = xc - hConvolvedSize; 
int ys = yc - hConvolvedSize;
const DImage &weight = stamp.weight;
for (int j=0; j< convolvedSize; ++j)
for (int i=0; i< convolvedSize; ++i)
  {
  double w_value = WorstImage(i+xs, j + ys);
  double res = R_conv_K(i,j) - w_value + BackValue(i + xs,j + ys); 
  double w = weight(i,j);
  chi2 += res*res*w;
  // the number of non zero weight pixels is already in Stamps
  }
stamp.chi2 = chi2;
 if(chi2<0) {
   cerr << " KernelFit::StampChi2 WARNING xc,yc,chi2 = " << xc << "," << yc << "," << chi2 << endl;  
 } 
return chi2;
}


bool KernelFit::Solve(const Mat &m, const Vect &b, const Point &Center)
{
  int oldprec = cout.precision();
  cout << setprecision(10);
  ios::fmtflags  old_flags = cout.flags(); 
  cout << resetiosflags(ios::fixed) ;
  cout << setiosflags(ios::scientific) ;
  if (solution.size()!=(unsigned int)mSize) solution.resize(mSize);
  bool could_solve;

#define NEWWAY
#ifdef NEWWAY
  double dchi2 = 0;
    /* operate on a copy, to preserve m & b, in case we subtract outlier
       stamps later. 
       m&b are "const" BTW.*/
  Mat mprime(m);
  Vect bprime(b);
  // we can use cholesky because m is (should be) posdef 
  could_solve = (cholesky_solve(mprime,bprime, "U") == 0);
  cout << " Kernel Inversion: " << could_solve << endl;
  if (!could_solve) return 0;
  if (optParams.KernVar.Nterms()>1 && optParams.UniformPhotomRatio) 
//  use Lagrange multipliers as a modification of the "free" solution
    {
      cholesky_invert(mprime,"U");
      int nc = optParams.KernVar.Nterms() -1; // number of constraints
      int nKern = Basis.size();
      cout <<" Integral of kernel is forced to be constant." << endl ;
      Mat C(nc, mSize); // the matrix of constraints
      for (int ik=0; ik < nKern; ++ik)
	{
	  double kern_int = Basis[ik].sum();
	  for (unsigned ic = 1; ic < optParams.KernVar.Nterms(); ++ic)
	    {
	      int ip = KernIndex(ik,ic);
	      int jp = ic-1;
	      C(jp,ip) = kern_int;
	    }
	}
      Mat mm1C(mprime * C);
      Mat D(C.transposed()*mm1C);
      Vect lambda(C.transposed()*bprime); // contains D*lambda for now (will be lambda once solved for)
      Vect Dlambda(lambda);
      if (cholesky_solve(D,lambda, "U") != 0)
	{
	  cout << " problems when solving for constraints in kernel fit.... " << endl;
	  return 0;
	}
      dchi2 = lambda*Dlambda;
      cout << "dchi2 for constant photom ratio =" << dchi2 << endl;
      bprime -= mm1C*lambda;
    }
    if (could_solve)
      {
	solution.resize(mSize);
	for(unsigned int i=0; i<mSize;i++) {
	  solution[i]=bprime(i);
	}
      }
#else /*OLDWAY */

  //DEBUG
 // i.e. the integral of the kernel (photometric ratio) is constant over the image

    if (optParams.UniformPhotomRatio) {
    // use Lagrange multipliers technique
    cout <<" Integral of kernel is assumed to be constant." << endl ;
    int nKern = Basis.size();
    int nc = optParams.KernVar.Nterms() -1; // number of constraints
    int totSize = mSize + nc;
    Mat mprime(totSize,totSize);
    Vect bprime(totSize);
    for (unsigned k=0; k<mSize; ++k) bprime(k) = b(k);
    for (unsigned i=0; i<mSize; ++i)
      for (unsigned j=0; j<mSize; ++j)
	{
	  mprime(i,j) = m(i,j);
	}
    for (int ik=0; ik < nKern; ++ik)
      {
	double kern_int = Basis[ik].sum();
	for (unsigned ic =1; ic < optParams.KernVar.Nterms(); ++ic)
	  {
	    int ip = KernIndex(ik,ic);
	    int jp = mSize+ic-1;
	    mprime(ip,jp) = mprime(jp,ip) = kern_int;
	  }
      }
    // DEBUG
    //mprime.writeFits("A.fits");

    /* have to use a lin eq. solver that accomodates non posdef
       matrices : mprime is NOT posdef. */
    could_solve = (general_solve(mprime, bprime,false /* no inverse */, "U") == 0);
    cout << " Kernel Inversion " << could_solve << endl;
    if (could_solve)
      {
	solution.resize(mSize);
	for(size_t i=0; i<mSize;i++)
	  solution[i]=bprime(i);
      }
    cout << " lambda ";
    for (int k = bprime.Size()-nc; k < bprime.Size(); ++k)
      cout << bprime(k) << endl;

  }// end of if (constant kernel integral) ...
  else {// NO CONSTRAINT on the variations of kernel integral
    /* operate on a copy, to preserve m & b, in case we subtract outlier
       stamps later. 
       m&b are "const" BTW.*/
    Mat mprime(m);
    Vect bprime(b);
    // we can use cholesky because m is (should be) posdef 
    could_solve = (cholesky_solve(mprime,bprime, "U") == 0);
    cout << " Kernel Inversion: " << could_solve << endl;
    if (could_solve)
      {
	solution.resize(mSize);
	for(unsigned int i=0; i<mSize;i++) {
	  solution[i]=bprime(i);
	}
      }
  }
#endif


  Kernel kernel_at_center( optParams.HKernelSize, optParams.HKernelSize);
  KernCompute(kernel_at_center, Center.x, Center.y);
  cout << " Kernel characteristics "; kernel_at_center.dump_info();

  if (optParams.BackVar.Nterms()) 
    {
      cout << setprecision(10);
      cout << " Differential background: " << endl;
      for (size_t ib=0; ib< optParams.BackVar.Nterms(); ++ib) 
	cout << solution[BackIndex(ib)] << " " ;
      cout << endl;
   }
  cout << setprecision(oldprec);
  cout.flags(old_flags);

  StoreScore("kfit","dchi2",dchi2); // dchi2 for constant integral

  return could_solve;
}

double KernelFit::SepBackValue(const double&x, const double &y) const
{
  double val=0.0;

  for (size_t ib=0; ib < optParams.SepBackVar.Nterms(); ++ib)
    {
      val += diffbackground[ib]*optParams.SepBackVar.Value(x,y,ib);
    }
  return val;
}


double KernelFit::BackValue(const double&x, const double &y) const
{
  double val=0.0;

  for (size_t ib=0; ib < optParams.BackVar.Nterms(); ++ib)
    {
      val += solution[BackIndex(ib)]*optParams.BackVar.Value(x,y,ib);
    }
  return val;
}


//! Carry out kernel fit.

int KernelFit::DoTheFit(ImagePair &ImPair)
{
  StarMatchList matchList;
  int nobj = MakeObjectList(ImPair, 1, matchList);
  cout << " We have " << nobj << " objects to fit with" << endl;

  if ( nobj == 0 ) {
    throw(PolokaException("No stars for psfmatch"));
  }
  if (nobj < 10) cerr << " WARNING: Less than 10 common objects for kernelfit "<< endl;
  cout << " Frame limits for the fit: " <<  ImPair.CommonFrame() << endl;

  double bestSeeing = ImPair.BestSeeing();
  double worstSeeing = ImPair.WorstSeeing();
  optParams.OptimizeSizes(bestSeeing, worstSeeing);
  cout << " Max stamps = " << optParams.MaxStamps << endl;
  optParams.OptimizeSpatialVariations(min(optParams.MaxStamps,int(matchList.size())));

 //cook up a plausible kernel to propagate weight map if bestimage
  Kernel guess(optParams.HKernelSize);
  if(worstSeeing>bestSeeing) {
    PolGaussKern(guess, sqrt(worstSeeing*worstSeeing - bestSeeing*bestSeeing), 0, 0);
    guess *= 1./guess.sum();
  }else{ // use delta
    SetDelta(guess);
  }

  double sexPhotomRatio = 1./MedianPhotomRatio(matchList);

  guess *=  sexPhotomRatio;
  StampList bestImageStamps(ImPair, matchList, optParams.HStampSize, guess, optParams.MaxStamps);

  double InitialSig2Noise = bestImageStamps.Sig2Noise();
  cout << " Total s2n of objects initially used for the fit " 
       << InitialSig2Noise << endl;
  StoreScore("kfit","is2n",InitialSig2Noise);

  nstamps = bestImageStamps.size();
  
  if (!nstamps)
    {
      cerr << " No Stars to fit the kernel ... : giving up " << endl;
      return 0;
    }

  cout << " starting kernel fit with 1/2 StampSize " << optParams.HStampSize 
       << " 1/2 KernelSize  " << optParams.HKernelSize 
       << " convolved Stamp size " << optParams.ConvolvedSize() << endl; 
  cout << " with " << bestImageStamps.size() << " stamps." << endl;
  clock_t tstart = clock();

  FitDifferentialBackground(ImPair, 3.0); // does nothing if datacards don't tell it explicitly

  // handle worst: it is the actual worst if diff. background is not fitted separatly,
  // we have to subtract diff. background if it was fitted:
  const Image *worstImage=NULL;
  auto_ptr<Image> del_worst(NULL); // to ensure deletion of temporary image, if needed.
  if (optParams.SepBackVar.Nterms() != 0)
    {
      Image *tmp = new Image(ImPair.WorstImage());
      optParams.SepBackVar.ApplyToImage(*tmp, -1, diffbackground);
      worstImage = tmp;
      del_worst.reset(tmp);
    }
  else
    {
      worstImage = &ImPair.WorstImage();
    }
  
  Mat basisTransfo;
  BasisFill(&basisTransfo);

/* allocate the convolved stamps space */
  FitWorkSpace mySpace(*this);
  work = &mySpace;

  Mat m(mSize,mSize);
  Vect b(mSize);

  // fill the LS problem matrix and RHS vector m and b (long)
  for (StampIterator si = bestImageStamps.begin(); 
       si != bestImageStamps.end(); ++si)
    OneStampMAndB(*si, *worstImage, m, b);

  if (0)
    { // DEBUG
    clock_t tstart = clock();
    Mat mnoise(mSize, mSize);
    for (StampIterator si = bestImageStamps.begin(); 
	 si != bestImageStamps.end(); ++si)
      OneStampNoiseMatrix(*si, *worstImage, mnoise);
    //    m -= mnoise;clock_t tend = clock();
    clock_t tend = clock();
    cout << "CPU for mnoise " <<  float(tend- tstart)/float(CLOCKS_PER_SEC) << endl;

    mnoise.writeFits("mnoise.fits");
  }
  

  cout << " finished computation of m and b" << endl;

  //now iterate stamp filtering and refitting until all stamps are acceptable
  int dropped;
  vector<double> chi2s(nstamps);
  int oldprec = cout.precision();
  cout << setprecision(10);
  // A REGARDER
  cout << " BackValue(0,0)" << BackValue(0,0) << endl;
  int npixtot=0;
  double chi2_tot=0;
  Point center = ImPair.CommonFrame().Center(); // for printouts
  do // iterations on stamp clipping.
    {
      if (!Solve(m,b, center))
	{
	  fitDone=false;
	  cerr << " KernelFit: Inversion failed  " << endl;
	  break;
	}
      fitDone = true;
      // compute stamps chi2 with this solution
      int i=0;
      npixtot = 0;
      chi2_tot = 0;
      for (StampIterator si = bestImageStamps.begin(); si != bestImageStamps.end(); ++si, ++i)
       {
	 Stamp &stamp = *si;
	 double chi2_stamp = StampChi2(stamp, *worstImage);
	 chi2s[i] = chi2_stamp/stamp.nActivePix; /* also assigns stamp.chi2 */
	 chi2_tot += chi2_stamp;
	 npixtot += stamp.nActivePix;
	 // cout << " xc yc chi2 " << stamp.xc << ' ' << stamp.yc << ' ' << chi2s[i] << endl;
	 // if (stamp.star) stamp.star->dump();
       }  
      double mean,sigma,median;
      mean_median_sigma(&(chi2s[0]),nstamps,mean,median,sigma);
      cout << " chi2/dof per stamp : mean median sigma " << mean << ' ' << median << ' ' << sigma << endl;
      cout << " chi2, ndof, chi2/ndof " << chi2_tot << ' ' << npixtot-mSize << ' ' << chi2_tot/(npixtot - mSize) << endl;
      dropped = 0;
     
      // this trick is used to fit the kernel using a predefined catalog
      // of objets to compute the kernel
      if (getenv("NOFILTERING")) break;
      for (StampIterator si = bestImageStamps.begin(); si != bestImageStamps.end(); )
	{
	  Stamp &stamp = *si; 
	  /* cut too large chi2's and negative ones (due to dead columns)*/
	  double chi2 = stamp.chi2/stamp.nActivePix;
	  if (chi2 > median + 4*sigma || chi2 < 0) 
	    {
	      cout << " delete : xc yc chi2 npix chi2/npix " << stamp.xc << ' ' << stamp.yc << ' ' 
		   << stamp.chi2 << ' ' << stamp.nActivePix << ' ' << chi2 << endl;
	      //if (stamp.star) stamp.star->dump();
	      // drop it :
	      dropped++;
	      Mat mStamp(mSize,mSize);
	      Vect bStamp(mSize);
	      OneStampMAndB(stamp, *worstImage,mStamp, bStamp);
	      m -= mStamp;
	      b -= bStamp;
	      si = bestImageStamps.erase(si);
	    }
	  else ++si;
	}
      if (dropped) cout << " dropped " << dropped << " stamps : refitting " << endl;
      nstamps -= dropped;
    } 
  while (dropped != 0);   // done with stamp clipping

  //DEBUG
  if (0)
    {
      Mat eigenVects(m);
      Vect eigenVals(b);
      //      symetric_diagonalize(eigenVects,eigenVals,"L");
      DiagonalizeRealSymmetricMatrix(m, eigenVects, eigenVals);
      cout << "lowest eigenvalues " << eigenVals(0) << ' ' << eigenVals(1) << ' ' << endl;
      StoreScore("kfit","lammin0",eigenVals(0));
      StoreScore("kfit","lammin1",eigenVals(1));
      StoreScore("kfit","lammin2",eigenVals(2));
    }

  if (!fitDone) return 0;
  if (0) {// STUDIES
    Mat groups;
    ParameterGroups(groups); // groups is expressed ignoring orthogonalisation...
    int nt = optParams.KernVar.Nterms(); 
    Mat paramTransfo(Basis.size()*nt, Basis.size()*nt);
    for (unsigned jk=0; jk<Basis.size(); ++jk)
      for (unsigned js=0; js<nt; ++js)
	for (unsigned ik=0; ik<Basis.size(); ++ik)
	  paramTransfo(nt*ik+js, nt*jk+js) = basisTransfo(ik,jk);
    Mat mt(paramTransfo.transposed()*m*paramTransfo);
    Vect bt(paramTransfo.transposed()*b);    
    vector<double> dchi2 = DeltaChi2(mt,bt,groups);
    StoreScore("kfit","dchi2s", dchi2);
  }

  // m.writeFits("A.fits");

  // mostly printouts from now on
  cout  << " finished the fit  with " << bestImageStamps.size() << " stamps " << endl;
  { // DEBUG
    char *fileName = getenv("DUMP_FIT_LIST");
      if (fileName)
	{
	  ofstream l(fileName);
	  bestImageStamps.front().star->WriteHeader_(l);
	  l << "# chi2 : \n# end " << endl;
	  for (StampIterator si = bestImageStamps.begin(); si != bestImageStamps.end(); ++si)
	    {
	      (*si).star->writen(l); 
	      l << ' ' << (*si).chi2 << endl;
	    }
	}
  }

  double FinalSig2Noise = bestImageStamps.Sig2Noise();
  cout << " Total s2n of objects eventually used for the fit " 
       << FinalSig2Noise << endl;
  StoreScore("kfit","fs2n",FinalSig2Noise);
  StoreScore("kfit","nparams",mSize);

  double mean,median,sigma;
  mean_median_sigma(&(chi2s[0]), nstamps, mean, median, sigma);
  chi2 = chi2_tot/(npixtot - mSize);
  cout << " final chi2 ndof chi2/dof " << chi2_tot << ' ' << npixtot - mSize << ' ' << chi2 << endl;
  StoreScore("kfit","kchi2",chi2_tot);
  StoreScore("kfit","kndof",npixtot - mSize);
  StoreScore("kfit","knstamps", nstamps);


  Kernel kernel_at_center( optParams.HKernelSize, optParams.HKernelSize);
  KernCompute(kernel_at_center, center.x, center.y);
  kernAtCenterSum = kernel_at_center.sum();
  double mxx, myy, mxy;
  kernel_at_center.moments(mxx,myy,mxy);
  StoreScore("kfit","kmxx",mxx);  
  StoreScore("kfit","kmyy",myy);
  StoreScore("kfit","kmxy",mxy);
  StoreScore("kfit","kmxy",mxy);
  StoreScore("kfit","kint",kernAtCenterSum);
  StoreScore("kfit","k2int", kernel_at_center.sum2());
  StoreScore("kfit","kvarn", optParams.KernVar.Nterms());
  StoreScore("kfit","hksize",optParams.HKernelSize);
  StoreScore("kfit","hssize",optParams.HStampSize);
  
  cout << setprecision(oldprec);

  // there are historical scripts that grep this line in the log :
  cout << "PsfMatch_SUMMARY_best_worst_kernelphotomratio_sextractorratio_nstamps_chi2/dof "
       << ImPair.Best()->Name() << ' ' << ImPair.Worst()->Name() << ' ' 
       << kernAtCenterSum << ' ' 
       << sexPhotomRatio << ' '
       << nstamps << ' ' 
       << Chi2() << ' ' 
       << endl;

  clock_t tend = clock();
  cout << "CPU for the kernel fit " <<  float(tend- tstart)/float(CLOCKS_PER_SEC) << endl;
  return 1;
}


/************************************   KernelFit I/O's *******************************/


void KernelFit::read(const string& FileName)
{
  ifstream s(FileName.c_str());
  if (!s)
    throw(PolokaException(" KernelFit::read: cannot read file "+FileName)); 
  read(s);
}


void KernelFit::read(istream& stream)
{
  string tmp_str;
  int version;
  try {
    stream.exceptions (ios::failbit); // allow exceptions to be raised
    stream >> tmp_str >> version;
    if (tmp_str != "[KernelFit]") throw(PolokaException(" Do not recognize a kernel in this file"));
    if (version <=1)
      { // these things are not used any longer
	double BestImageBack, WorstImageBack, SkyVarianceWorstImage, 
	  WorstImageGain;
	read_member(stream, tmp_str, BestImageBack);
	read_member(stream, tmp_str, WorstImageBack);
	read_member(stream, tmp_str, SkyVarianceWorstImage);
	read_member(stream, tmp_str, WorstImageGain);
      }
    read_member(stream, tmp_str, kernAtCenterSum);
    read_member(stream, tmp_str, optParams);
    read_member(stream, tmp_str, mSize);
    read_member(stream, tmp_str, solution);
    read_member(stream, tmp_str, diffbackground);
    read_member(stream, tmp_str, chi2);
    read_member(stream, tmp_str, nstamps);
  }
  catch (const std::ios_base::failure& error)
    {
      throw(PolokaException(error.what()));
    }
  BasisFill(NULL);  // behavior depends on optParams, retrieved from the file.
  fitDone = true;
}

void KernelFit::write(const string &FileName) const
{
  if (!fitDone) 
    {
      cout << "ERROR: KernelFit::write : we don't write a Kernel that was not fitted (yet) " << endl; 
      return;
    }    
  ofstream s(FileName.c_str());
  write(s);
  s.close();
}


void KernelFit::write(ostream& stream) const
{
  /* versions : 
     0 = "old unweighted code"
     1 = new weighted code (april 2007). same format as version 0
     2 = after cleanup (march 2010) : some cores disapeared
  */
  stream << "[KernelFit] " << 2;
  stream << setprecision(18);
  /* these things ahve disappeared :
    write_member(stream, "BestImageBack", BestImageBack);
    write_member(stream, "WorstImageBack", WorstImageBack);
    write_member(stream, "SkyVarianceWorstImage", SkyVarianceWorstImage);
    write_member(stream, "WorstImageGain", WorstImageGain);
  */
  write_member(stream, "KernAtCenterSum", kernAtCenterSum);
  write_member(stream, "optParams", optParams);
  write_member(stream, "mSize", mSize);
  write_member(stream, "solution", solution);
  write_member(stream, "diffbackground", diffbackground);
  write_member(stream, "chi2", chi2);
  write_member(stream, "nstamps", nstamps);
}


#ifdef USE_ROOT
/*
RUN_ROOTCINT
LINKDEF_CONTENT : #pragma link C++ class OptParams+;
LINKDEF_CONTENT : #pragma link C++ class XYPower+;
LINKDEF_CONTENT : #pragma link C++ class KernelFit+;

*/
#include "root_dict/kernelfitdict.cc"

#endif /* USE ROOT */
