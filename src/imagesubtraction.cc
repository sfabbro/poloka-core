#include <algorithm>
#include <fstream>
#include <cassert>

#include "imagesubtraction.h"
#include "kernelfit.h"
#include "fileutils.h"
#include "fitsimage.h"
#include "imageback.h"
#include "gtransfo.h"
#include "apersestar.h"
#include "toadscards.h"
#include "datacards.h"
#include "quali_box.h"
#include "convolution.h"
#include "toadscards.h"
#include "detection.h"
#include "allreducedimage.h"
#include "imagepair.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

string SubtractedName(const string &RefName, const string &NewName)
{
  return "Sub_" + RefName + '_' + NewName;
}

ImageSubtraction::ImageSubtraction(const string &Name,  
				   const ReducedImageRef RefImage,
				   const ReducedImageRef NewImage,
				   const bool NoSwap)
  : ReducedImage(Name),
    KernelFitter(RefImage, NewImage, NoSwap),
    Ref(RefImage), New(NewImage)
{
  if (!FileExists(Name))
    Create("here");
}

ImageSubtraction::ImageSubtraction(const string &Name,
				   const ReducedImageRef RefImage,
				   const ReducedImageRef NewImage, 
				   const KernelFitter &APreviousFit)
  : ReducedImage(Name),
    KernelFitter(APreviousFit),
    Ref(RefImage), New(NewImage)
{
  Create("here");
}

ImageSubtraction::ImageSubtraction(const string &Name)
  : ReducedImage(Name),
    KernelFitter(KernelName())
{
  vector<string> dbNames;
  DecomposeString(dbNames, Name, "_");
  if (dbNames.size() == 3)
    {
      Ref = ReducedImageNew(dbNames[1]);
      New = ReducedImageNew(dbNames[2]);
      cout << " ImageSubtraction: ref " << Ref->Name() << " " << New->Name() << endl;
    }
  else
    {
      cerr << " ImageSubtraction: " << Name << " is an old subtraction, no persistence\n";
    }
}

ReducedImageRef ImageSubtraction::Clone() const
{
  return new ImageSubtraction(*this);
}

static double sqr(double x) { return x*x; }

// this routine creates a kind of ntuple that allows to study 
// the subtraction quality  (the famous Delphine's plot)
static void quali_plots(const ReducedImageRef Ref,
			const ReducedImageRef New, 
			const FitsImage *SubImage,
			const string& LogName)
{
  const Image *imArray[3];
  FitsImage refImage(Ref->FitsName());
  FitsImage newImage(New->FitsName());
  imArray[0] = &refImage;
  imArray[1] = &newImage;
  imArray[2] = SubImage;
  // the last one is zero because the kernelfit is supposed to achieve that
  double fond[3] = {Ref->BackLevel(), New->BackLevel(), 0};
  // slightly wrong sigma
  double sigmaFond[3] = {Ref->SigmaBack(), New->SigmaBack(),
			 sqrt(sqr(Ref->SigmaBack())+sqr(New->SigmaBack()))} ;
  
  SEStarList refList(Ref->CatalogName());
  SEStarList newList(New->CatalogName());
  ofstream qual_log(LogName.c_str());
  GtransfoIdentity ident;
  Test_Qual_N_same(3, *SubImage, refList, newList,
		   &ident, imArray, fond, sigmaFond, 10, qual_log);
}                 

bool ImageSubtraction::MakeFits()
{
  string fileName = FitsName();
  if (FileExists(fileName)) return true;
  cout << " ImageSubtraction: making " << fileName << endl;
  if (!Ref->MakeFits() || !New->MakeFits())
    {
      cerr << " ImageSubtraction: " << Name()
	   << " could not make images of (both) subtraction terms :"
	   << Ref->Name() << " and " << New->Name() << endl;
      return false;
    }

  // fit the convolution kernel
  if (solution.empty()) DoTheFit();

  // build the subtracted image
  FitsHeader *hnew = new FitsHeader(New->FitsName());
  FitsImage subFits(fileName, *hnew);
  delete hnew;
  Image &subImage = subFits;

  // by convention the subtraction photometric scale is the same as the ref
  //  - photomRatio is New/Ref
  //  - Important point : the same trick has to be applied to the variance computation.
  // convention changed: now subtraction is always on new (photomratio is still new/ref)
  double photomRatio = PhotomRatio();
  if (RefIsBest())
    {
      cout << " ImageSubtraction: " << New->Name() 
	   << " - Kernel*" << Ref->Name() << endl;
      {
	FitsImage worstFits(New->FitsName());
	subImage = worstFits;
      }
      {
	FitsImage bestFits(Ref->FitsName());
	Image bestImageConv(bestFits.Nx(), bestFits.Ny());
	ImageConvolve(bestFits, bestImageConv);
	AddBackground(bestImageConv);
	subImage -= bestImageConv;
      }
      // "*" is far faster than "/"
      //Pixel factor = 1./photomRatio;
      //subImage *= factor;
    }
  else 
    {
      cout << " ImageSubtraction: " << New->Name() 
	   << "*Kernel - " << Ref->Name() << endl;
      {
	FitsImage bestFits(New->FitsName());
	ImageConvolve(bestFits, subImage);
	AddBackground(subImage);
      }
      {
	FitsImage worstFits(Ref->FitsName());
	subImage -= worstFits;
      }
      // "*" is far faster than "/"
      Pixel factor = 1./photomRatio;
      subImage *= factor;
    }
  // force floating point values (not sure it is useful)
  subFits.AddOrModKey("KERNREF", Best()->Name(), " name of the seeing reference image");
  subFits.AddOrModKey("SUBREF", Ref->Name(), " name of the ref image");
  subFits.AddOrModKey("SUBNEW", New->Name(), " name of the new image");
  subFits.AddOrModKey("KERNCHI2", Chi2(), " chi2/dof of the kernel fit");
  subFits.AddOrModKey("PHORATIO", photomRatio, "flux = PHORATIO * flux(KERNREF)");
  double mjd= New->ModifiedJulianDate();
  subFits.AddOrModKey("MJD-OBS", mjd,"Modified Julian Date");

  // subtract an ImageBack
  FitsImage *pweight = NULL;  
  if (MakeWeight())
    {
      pweight = new FitsImage(FitsWeightName());
      cout << " ImageSubtraction: using weights : " << FitsWeightName() << endl;
    }
  
  string datacards = DefaultDatacards();
  int backMesh = 64;
  if (FileExists(datacards))
    {
      DataCards cards(datacards);
      if (cards.HasKey("KFIT_MESHSTEP"))
	{
	  backMesh = cards.IParam("KFIT_MESHSTEP");
	}
    }


  ImageBack b(subImage, backMesh, pweight);  
  Image *back = b.BackgroundImage();
  subImage -= *back ;
  delete back;
  subFits.AddOrModKey("BACK_SUB", true, 
		      pweight? " subtracted weighed background" :
		      " subtracted unweighted background");
  subFits.AddOrModKey("BACKMESH", backMesh,
		      " mesh size used for back computation ");

  if (pweight)
    {
      pweight->Simplify(1e-30);
      subImage *= *pweight;
      delete pweight;
    }

  //subFits.SetWriteAsFloat();
  subFits.ModKey("BITPIX",16);

  if (FileExists(datacards))
    {
      DataCards cards(datacards);
      if (cards.HasKey("SUB_QUALI_TUPLE"))
	{
	  string qualfile = "qual.list";
	  qualfile = cards.SParam("SUB_QUALI_TUPLE");
	  quali_plots(Ref, New, &subFits, AddSlash(Dir()) + qualfile);
	}
    }
  

  // set some score values
  SetSeeing(LargestSeeing());
  SetBackLevel(0.0); // by construction
  SetSaturation(Worst()->Saturation()," worst image saturation");
  double factor=1;
  if (!RefIsBest()) factor=1/photomRatio;
  SetSigmaBack(sqrt(sqr(Best()->SigmaBack()*photomRatio) + sqr(Worst()->SigmaBack()))*factor);
  SetReadoutNoise(sqrt(sqr(photomRatio*Best()->ReadoutNoise()) + sqr(Worst()->ReadoutNoise()))*factor);
  SetUsablePart(CommonFrame());
  SetOriginalSkyLevel((photomRatio*Best()->OriginalSkyLevel() + Worst()->OriginalSkyLevel())*factor
		      ," sum of sky levels of the subtraction terms");

  //double zero_point = Ref->AnyZeroPoint();
  //cout << " ImageSubtraction: zero point taken from reference: " << zero_point << endl;
  //SetZZZeroP(zero_point, "as taken from reference");
  
  return true;
}


static void add_model_noise(const SEStar *Star, 
			    const double& Gain,
			    Image &Variance)
{
 double mxx, myy, mxy;
  const AperSEStar* ap = dynamic_cast<const AperSEStar*>(Star);
  if (ap)
    {
      mxx = ap->gmxx;
      myy = ap->gmyy;
      mxy = ap->gmxy;
    }
  else 
    {
      mxx = Star->Mxx();
      myy = Star->Myy();
      mxy = Star->Mxy();
    }

  // account for object noise using a gaussian "model"
  double det = mxx * myy - sqr(mxy);
  double wxx = myy/det;
  double wyy = mxx/det;
  double wxy = mxy/det;
  double factor = Star->flux/Gain * (2*M_PI*sqrt(det));
  double hwidth = 3*Star->A() + 2;
  int xstart = max(0, int(floor(Star->x - hwidth)));
  int ystart = max(0, int(floor(Star->y - hwidth)));
  int xend = min(int(ceil(Star->x + hwidth)), Variance.Nx());
  int yend = min(int(ceil(Star->y + hwidth)), Variance.Ny());

  for (int j=ystart; j<yend; ++j)
    {
      double y = j - Star->y;
      for (int i=xstart; i<xend; ++i)
	{
	  double x = i - Star->x;
	  Variance(i,j) += factor * exp(-0.5*(wxx*x*x + wyy*y*y + 2*wxy*x*y));
	}
    }
}

// this routine computes the weight map of the subtraction:
//   - weights from best image
//   - account for satur in this image
//   - go to variances
//   - convolve (by the squared kernel form the image kernel fit)
//   - take weights from worst and account for satur in worst.
//   - combine both maps to get total variance
//   - normalize with the same convention as the image
//   - go back to weights
//   - zero the sides of the weight map
bool ImageSubtraction::MakeWeight()
{
  string fileName = FitsWeightName();
  if (FileExists(fileName)) return true;
  if (solution.empty()) DoTheFit();

  cout << " ImageSubtraction: making weight image  " << fileName << endl;
  string datacards = DefaultDatacards();
  int subVarianceType = 0;
  int subMaskDilate = 0;
  double subMaskSatur = 1;
    
  if (FileExists(datacards))
    {
      DataCards cards(datacards);
      if (cards.HasKey("SUB_VARIANCE_TYPE"))
	subVarianceType = cards.IParam("SUB_VARIANCE_TYPE");
      if (cards.HasKey("SUB_MASK_DILATE"))
	subMaskDilate = cards.IParam("SUB_MASK_DILATE");
      if (cards.HasKey("SUB_MASK_SATUR"))
	subMaskSatur = cards.DParam("SUB_MASK_SATUR");
    }
    
  FitsHeader refHead(Ref->FitsName());
  FitsImage subWeightFits(FitsWeightName(), refHead);
  CommonFrame().WriteInHeader(subWeightFits);
  Image& subWeightImage = subWeightFits;

  { // to open temporarily Best weight image.
    FitsImage bestWeightFits(Best()->FitsWeightName());
    subWeightImage = bestWeightFits;
  }
  
  // weight saturated pixels in Best to zero
  if (Best()->HasSatur())
    {
      FitsImage bestSaturFits(Best()->FitsSaturName());
      cout << " ImageSubtraction: masking from " << Best()->Name() << " saturation map\n";
      subWeightImage *= (1 - bestSaturFits);
    }
  else
    {
      cerr << " ImageSubtraction: no saturation masking from " 
	   << Best()->Name() << endl;
    }

  const Pixel *pend = subWeightImage.end();

  // weight down to zero pixels close to satur on best image
  if (subMaskSatur < 1 && RefIsBest())
    {
      FitsImage bestFits(Best()->FitsName());
      cout << " ImageSubtraction: Masking pixels at " 
	   << subMaskSatur*100 << "% of saturation from " << Best()->Name() << endl;
      Pixel *pim = bestFits.begin();
      Pixel satur = Best()->Saturation() * subMaskSatur;
      for (Pixel *pw=subWeightImage.begin() ; pw < pend ; ++pim, ++pw)
	if (*pim >= satur) *pw = 0;
    }
  
  // add a small constant to weights so that variances of 
  // zero weight pixels remain finite and go to variances
  Pixel eps = subWeightImage.MaxValue() * 1e-10;
  {
    for (Pixel *pw = subWeightImage.begin(); pw < pend; ++pw)
      *pw = 1. / (*pw + eps);
  }

  // Variance type for subtraction:
  //   0: no extra variance var(sub) = sky(best) * k^2 + sky(worst) default
  //   1: var(sub) = k^2 * (sky(best) + best) + sky(worst) + worst
  //      includes Poisson noise from objects as "data": bias photometry but useful when many bright stars
  //   2: var(sub) = k^2 * (sky(best) + model(best)) + sky(worst) + model(worst)
  //      includes poisson noise from objects as Gaussian model. Might be less biased than 
  //	  previous, but usually very wrong for crowded fields and sextractor based catalogs

  Pixel factfmax = 10.;
  if (subVarianceType == 1)
    {
      cout << " ImageSubtraction: weighting with data as Poisson noise\n";
      Pixel fact = 1. / Best()->Gain();
      Pixel back = Best()->BackLevel();
      FitsImage bestFits(Best()->FitsName());
      for (Pixel *pim=bestFits.begin(), *pw=subWeightImage.begin(); pw < pend; ++pim, ++pw)
	*pw += (*pim - back) * fact;
    } 
  else if (subVarianceType == 2)
    {
      cout << " ImageSubtraction: weighting with Gaussian models as Poisson noise\n";
      Pixel gain = Best()->Gain();
      Pixel minfmax = factfmax * sqr(Best()->SigmaBack());
      SEStarList bestList;
      if (Best()->HasAperCatalog()) bestList.read(Best()->AperCatalogName());
      else bestList.read(Best()->CatalogName());
      for (SEStarCIterator it=bestList.begin(); it != bestList.end(); ++it)
	if ((*it)->Fluxmax() > minfmax) add_model_noise(*it, gain, subWeightImage);
    } 
  else
    {
      cout << " ImageSubtraction: weight is built from sky noise\n";
    }

  // convolve with squared image kernel
  {
    Image toConvolve(subWeightImage);
    VarianceConvolve(toConvolve, subWeightImage);
  }

  // load the other weight map
  FitsImage worstWeightFits(Worst()->FitsWeightName());
  Image& worstWeightImage = worstWeightFits;

  // zero saturated pixels in Worst weight map
  if (Worst()->MakeSatur())
    {
      FitsImage worstSaturFits(Worst()->FitsSaturName());
      cout << " ImageSubtraction: masking from " << Worst()->Name() << " saturation map\n";
      worstWeightImage *= (1 - worstSaturFits);
    }
  else
    {
      cerr << " ImageSubtraction: no saturation masking from " 
	   << Worst()->Name() << endl;
    }

  if (subMaskSatur > 0 && !RefIsBest())
    {
      FitsImage worstFits(Worst()->FitsName());
      cout << " ImageSubtraction: masking pixels with satur > "
	   << subMaskSatur << " on " << Worst()->Name() << endl;
      Pixel *pim = worstFits.begin();
      Pixel *pww = worstWeightImage.begin();
      Pixel back = Worst()->BackLevel();
      for (Pixel *pws=subWeightImage.begin() ; pws < pend ; ++pim, ++pws, ++pww)
	if ((*pim-back) * sqrt(*pww) >= subMaskSatur) *pws = 0;
    }

  // add a small constant to weights so that variances of 
  // worse weights zero weight pixels remain finite and
  // now go to variances
  {
    eps = worstWeightImage.MaxValue() * 1e-10;
    const Pixel *pwend = worstWeightImage.end();    
    for (Pixel *pw = worstWeightImage.begin(); pw < pwend; ++pw) 
      *pw = 1. / (*pw + eps);

  }

  // add the poisson noise from objects from worst
  if (subVarianceType == 1)
    {
      Pixel fact = 1. / Worst()->Gain();
      Pixel back = Worst()->BackLevel();
      FitsImage worstFits(Worst()->FitsName());
      Pixel *pw = worstWeightImage.begin();
      const Pixel *pwend = worstWeightImage.end();     
      for (Pixel *pim=worstFits.begin(); pw < pwend; ++pw, ++pim)
	*pw += (*pim - back) * fact;
    }
  else if (subVarianceType == 2)
    {
      Pixel gain = Worst()->Gain();
      Pixel minfmax = factfmax * sqr(Worst()->SigmaBack());
      SEStarList worstList;
      if (Worst()->HasAperCatalog()) worstList.read(Worst()->AperCatalogName());
      else worstList.read(Worst()->CatalogName());

      for (SEStarCIterator it=worstList.begin(); it != worstList.end(); ++it)
	if ((*it)->Fluxmax() > minfmax)
	  add_model_noise(*it, gain, worstWeightImage);
    }

  // add best and worst weights, make a mask and dilate it
  {
    Pixel threshold = eps * 100;
    Image mask(subWeightImage.Nx(), subWeightImage.Ny());
    Pixel *pww = worstWeightImage.begin();
    Pixel *pm = mask.begin();
    for (Pixel *psw=subWeightImage.begin();  psw<pend; ++psw, ++pww, ++pm)
      {
	*psw = *pww / ((*psw * *pww) + 1);
	if (*psw < threshold) { *psw = 0; *pm = 1; }
	else *pm = 0;
      }
    if (subMaskDilate > 0)
      {
	cout << " ImageSubtraction: dilating bad pixels with " << subMaskDilate << " pixels\n";
	dilate_binary_image(mask, subMaskDilate);
	subWeightImage *= (1 - mask);
      }
  }

  // account for the fact that the subtraction is by convention
  // expressed in units of the ref. The weight map we have is
  // correctly normalized only if New was convolved (and
  // photometrically matched to ref by the kernel fit).  Apply the same
  // factor as the one applied to the image
  // convention has changed: units of New, so image unchanged
  if (!RefIsBest())
    {
      double ratio2 = sqr(PhotomRatio());
      cout << " ImageSubtraction: multiplying weight with with " << ratio2 << endl;
      subWeightImage *= ratio2;
    }

  // set to 0 the "side bands" where the convolution did not go
  Frame aframe(CommonFrame()); // whole image
  Kernel kern; // get a kernel to grab its size.
  KernAllocateAndCompute(kern, subWeightImage.Nx()/2., subWeightImage.Ny()/2);
  // size of the dead band due to variance convolution
  int bandx = kern.Nx()/2 -1;
  int bandy = kern.Ny()/2 -1;
  // shrink the frame
  aframe.CutMargin(bandx, bandy); 
  cout << " ImageSubtraction: masking frame for weights " << aframe;
  subWeightImage.Masking(aframe, 0.);
  float wmin, wmax;
  subWeightImage.MinMaxValue(&wmin,&wmax);
  if (fabs(wmin - wmax) < 1e-10) 
    cerr << " ImageSubtraction: " << Name() << " Warning: weight seems to be zero everywhere\n";
  subWeightFits.PreserveZeros(); // 0 remain 0 on R/W operations
  subWeightFits.ModKey("BITPIX", 16); // 16 bits are enough

  return true;
}

bool ImageSubtraction::MakeDead()
{
  if (FileExists(FitsDeadName())) return true;
  ReducedImageList imList(false);
  imList.push_back(Ref);
  imList.push_back(New);
  bool return_value = BoolImageOr(imList, 
				  &ReducedImage::FitsDeadName, 
				  &ReducedImage::MakeDead, FitsDeadName());
  FitsImage deadFits(FitsDeadName(), RW);
  cout << " ImageSubtraction: frame for dead of " << Name() << " is " << UsablePart();
  double presum = deadFits.SumPixels();
  deadFits.Masking(UsablePart(),1);
  cout << " ImageSubtraction: ratio of masked dead pixels after masking "
       << deadFits.SumPixels()/presum << endl;
  return return_value;
}

bool ImageSubtraction::MakeSatur()
{
  if (FileExists(FitsSaturName())) return true;
  ReducedImageList imList(false);
  imList.push_back(Ref);
  imList.push_back(New);
  return BoolImageOr(imList,
		     &ReducedImage::FitsSaturName, 
		     &ReducedImage::MakeSatur, FitsSaturName());
}

bool ImageSubtraction::MaskSatur()
{
   if (!FileExists(FitsName()) || !FileExists(FitsSaturName())) return false; 
   FitsImage subFits(FitsName(),RW);
   FitsImage subSaturFits(FitsSaturName());
   subFits *= (1.-subSaturFits);
   return true;
 }

bool ImageSubtraction::MaskNullWeight()
{
 if (!FileExists(FitsName()) || !FileExists(FitsWeightName())) return false;
 FitsImage subFits(FitsName(),RW);
 FitsImage subWeightFits(FitsWeightName());
 subWeightFits.Simplify(1e-30);
 subFits *= subWeightFits;
 return true;
}


bool ImageSubtraction::MakeCosmic()
{
  if (FileExists(FitsCosmicName())) return true;
  MakeWeight();

  FitsImage subFits(FitsName());
  if (!subFits.IsValid())
    {
      cerr << " ImageSubtraction: missing image to make cosmic image\n";
      return false;
    }

  FitsImage subWeightFits(FitsWeightName(), RW);
  if (!subWeightFits.IsValid())
    {
      cerr << " ImageSubtraction: missing weight to count cosmics\n";
      return false;
    }
  subFits *= subWeightFits;
  Pixel mean, sigma;
  subFits.SkyLevel(&mean, &sigma);
  FitsImage subCosmicFits(FitsCosmicName(), (FitsHeader&)subFits);;
  subFits.Cosmics(sigma, mean, Seeing(), subCosmicFits);
  subCosmicFits.AddOrModKey("BITPIX",8);
  // update weights
  subWeightFits *= (1.- subCosmicFits);
  subWeightFits.AddOrModKey("COSMPIXS", true, "This weight accounts for (identified) cosmics");
  return true;
}

bool ImageSubtraction::MakeCatalog()
{
  if (HasCatalog()) return true;
  cout << " ImageSubtraction: detecting point sources on " << Name() << endl;
  DetectionList detections;
  return (ImageDetect(*this, detections, Ref) && 
	  detections.write(CatalogName()) == 1);
}
