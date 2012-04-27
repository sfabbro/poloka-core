#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#include "fitsimage.h"
#include "reducedimage.h"
#include "transformedimage.h"
#include "allreducedimage.h"
#include "fileutils.h"
#include "vutils.h"
#include "imagesum.h"
#include "reducedutils.h"
#include "photoratio.h"

static double sqr(const double &a) { return a*a; }

static bool IsAlignedWith(const ReducedImage& One, const ReducedImage& Two) {
  FitsHeader h1(One.FitsName());
  if (h1.HasKey("GEOREF")) {
    string r1 = h1.KeyVal("GEOREF");
    if (r1 == Two.Name())
      return true;
    FitsHeader h2(Two.FitsName());
    if (h2.HasKey("GEOREF")) {
      string r2 = h2.KeyVal("GEOREF");
      return (r1 == r2);
    }
  }
  return false;
}

static bool AreAligned(const ReducedImage& One, const ReducedImage& Two) {
  return IsAlignedWith(One,Two) || IsAlignedWith(Two,One);
}

Component::Component(ReducedImage *RI, const double &PhotomRatio, 
		     const WeightingMethod weightingMethod) 
  :  backVar(0.), seeing(2.5), globalWeight(1.), averageWeight(1.), Ri(RI)
{
  backVar = sqr(Ri->SigmaBack());
  if (backVar == 0) // sigmaBackRed missing :  compute it the hard way
    {
      cout << " Component: did not find Sigma back in image " << Ri->Name() 
	   << ", compute it the hard way " << endl;
      FitsImage image(Ri->FitsName());
      float backRed, sigBackRed;
      image.SkyLevel(&backRed,&sigBackRed);
      backVar = sqr(sigBackRed);
    }  

  seeing = Ri->Seeing();
  photomRatio = PhotomRatio; /* of this one w.r.t the same reference for 
				all involved images.*/
  SetGlobalWeights(weightingMethod);
}

void  Component::SetGlobalWeights(const WeightingMethod weightingMethod)
{
  globalWeight = 1.;
  switch (weightingMethod)
    {
      /* the global weights are proportional to (signal/noise)^2
	 However we should not account here for sky background because 
	 this is already in the local weights. The photometric ratios 
	 are also not included in the global weight, but treated separately
	 (this is to be documented in detail..) */
      
    case PointSourceOptimal :
      /* noise for a faint point source scales as (seeing*sigback). 
	 signal scales as photomRatio. ignore sigback (see above). */
      globalWeight = 1./sqr(seeing);
      break;
    case ExtendedSourceOptimal :
      // noise scales as sigback. Ignore sigback (see above)
      break;
    case NoGlobalWeighting :
      photomRatio = 1.;
      break;
    case NoWeightsAtAll :
      photomRatio = 1.;
      break;
    default:
      cerr << " Component: the provided WeightingMethod :" << weightingMethod 
	   << " is unknown " << endl;
      globalWeight = 0.;
      break;
    }
}


struct PixVal {
  double pixVal; // p
  double fluxVal; // f
  double localWeight; // w 
  double globalWeight; // w'
  double totalWeight; // w*w'
  double pixVar; // Var(p). Sometimes, it is != to 1/w (stacking options)
  double photomRatio; // phi
  bool keep;
};


void Component::dump(ostream &s) const
{
  // have to restore Formats from C
  int oldprec = s.precision();
  s << setprecision(6);
  s << "Component from " << Ri->Name() << " : " 
    << "backVar =" << backVar 
    << " seeing = " << seeing 
    << " photomratio = " << photomRatio 
    << " globalweight = " << globalWeight 
    << " averageweight = " << averageWeight 
    << endl;
  s << setprecision(oldprec);
}


/****************  Stacking datacards handling **********/

#include "datacards.h"
#include "toadscards.h"

struct DatStack {

  WeightingMethod weightingMethod;
  StackingMethod stackingMethod;
  PhotoRatioMethod scalingMethod;

  DatStack(const string &DatacardsFileName);
  void Print(ostream& s=cout)const;
};

const string WEIGHTING_METHOD = "WEIGHTING_METHOD";
const string STACKING_METHOD = "STACKING_METHOD";
const string PHOTOSCALING_METHOD = "PHOTOSCALING_METHOD";
const string BACKGROUND_METHOD = "BACKGROUND_METHOD";

DatStack::DatStack(const string &DatacardsFileName)
{
  weightingMethod = ExtendedSourceOptimal;
  stackingMethod = WeightedAverage;
  scalingMethod = TotalLeastSquares;

  if (FileExists(DatacardsFileName))
    {
      DataCards cards(DatacardsFileName);

      if (cards.HasKey(WEIGHTING_METHOD))
	weightingMethod = (WeightingMethod) cards.IParam(WEIGHTING_METHOD);
      else 
	cerr << " DatStack: using old fashioned datacards : missing " 
	     << WEIGHTING_METHOD;
      
      if (cards.HasKey(STACKING_METHOD))
	stackingMethod = (StackingMethod) cards.IParam(STACKING_METHOD);
      else
	cerr << " DatStack: using old fashioned datacards : missing " 
	     << STACKING_METHOD;

      if (cards.HasKey(PHOTOSCALING_METHOD))
	scalingMethod = (PhotoRatioMethod) cards.IParam(PHOTOSCALING_METHOD);
      else
	cerr << " DatStack: using old fashioned datacards : missing " 
	     << PHOTOSCALING_METHOD;

    }
  else
    {
      cerr << " DatStack: cannot read " << DatacardsFileName << endl;
    }
}


void DatStack::Print(ostream& s) const
{
  s << "weighting method : " << weightingMethod << endl
    << "stacking method : " << stackingMethod << endl
    << "photoscaling method : " << scalingMethod << endl;
}


/***************** ImageSum ********************/

ImageSum::ImageSum(const string &AName, ReducedImageList &Images,
		   const ReducedImage *PhotomReference,
		   const WeightingMethod AWMethod, 
		   const StackingMethod ASMethod,
		   const PhotoRatioMethod APMethod)
 : ReducedImage(AName)
{

  DatStack datstack(DefaultDatacards());

  if (AWMethod == WUnSet)
    weightingMethod = datstack.weightingMethod;
  else
    weightingMethod = AWMethod;

  if (ASMethod == SUnSet)
    stackingMethod = datstack.stackingMethod;
  else
    stackingMethod = ASMethod;

  if (APMethod == PUnSet)
    scalingMethod = datstack.scalingMethod;
  else
    scalingMethod = APMethod;

  cout << " ImageSum: constructing " << Name() << endl
       << " ImageSum: " << Images.size() << " of images provided\n";

  if (Images.empty()) return;

  double total_weight = 0;

  const ReducedImage *photomReference = PhotomReference;
  if (!photomReference) photomReference = Images.front();
  photomReferenceName = photomReference->Name();
  zero_point_ref = photomReference->AnyZeroPoint();
  cout << " ImageSum: " << photomReferenceName << " is photometric reference\n";

  for (ReducedImageIterator i = Images.begin(); i!= Images.end(); ++i)
    {
      ReducedImage *ri = *i;
      bool ok = true ;
      if (!ri->HasImage())
	{
	  cout << " ImageSum: building missing " << ri->FitsName() << endl;
	  ok =  ri->MakeFits();
	}
      if (!ok) continue ;
      double phRatio = 1;
      if (weightingMethod != NoGlobalWeighting && 
	  weightingMethod != NoWeightsAtAll) {
	double err;
	phRatio = PhotoRatio(*ri, *photomReference, err, 0, scalingMethod);
      }
      if (!ok) continue;
      Component current(ri, phRatio, weightingMethod);
      current.dump();
      components.push_back(current);

      total_weight += components.back().globalWeight;
    }

  if (components.empty())
    {
      cerr << " ImageSum: no components to build ImageSum " << Name() << endl
	   << " ImageSum: giving up\n";
      return;
    }      
  size_t nComponents = components.size();
  total_weight /= double(nComponents);
  /* normalize weights. This is just esthetical because the 
     result should not depend on a global factor applied to weights */
  if (total_weight == 0) 
    {
      cerr << " ImageSum: we have images to build "<< Name() 
	   << ", but they all have a null weight ... giving up\n";
      return;
    }
  for (size_t k=0; k<nComponents; ++k) 
    {
      components[k].globalWeight /= total_weight;
      //components[k].dump();
    }

  if (components.empty())
    cerr << " ImageSum: constructed ImageSum ("<< Name() 
	 << ") without anything to sum " << endl;
  //  cerr << " Now creating Dbimage " << endl;
  Create("here");
}

ReducedImageRef ImageSum::Clone() const
{
  return new ImageSum(*this);
}

// this is the constructor for an already existing image
ImageSum::ImageSum(const string &Name) : ReducedImage(Name)
{
  // should read here the saved file.
}

void ImageSum::dump(ostream &s) const
{
  s << " ImageSum \"" << Name() << "\" : made from ";
  for (ComponentCIterator i=components.begin(); i != components.end(); ++i)  
    s << (*i).Ri->Name() << ' ';
  s << endl;
}

string name_of_stackingMethod(const StackingMethod stackingMethod)
{
  switch (stackingMethod)
    {
    case WeightedAverage :
      return "WeightedAverage";
      break;

    case ClippedWeightedAverage :
      return "ClippedWeightedAverage";
      break;

    case Median : 
      return "Median";
      break;

    case AdaptiveWeightedAverage : 
      return "AdaptiveWeightedAverage";
      break;

    default :
      cerr << " stacking method: " << stackingMethod 
	   << " is out of bounds " << endl;
    }
  return " ";
}

string name_of_weightingMethod(const WeightingMethod weightingMethod)
{
  switch (weightingMethod)
    {
    case PointSourceOptimal :
      return "PointSourceOptimal";
      break;

    case ExtendedSourceOptimal :
      return "ExtendedSourceOptimal";
      break;

    case NoGlobalWeighting : 
      return "NoGlobalWeighting";
      break;

    case NoWeightsAtAll :
      return "NoWeightsAtAll";
      break;

    default :
      cerr << " weighting method " << weightingMethod 
	   << " is out of bounds " << endl;
    }
  return " ";
}


#include "vutils.h"
// may go in vutils

static bool PixVal_fluxVal_compare(const PixVal &one, const PixVal &two) 
{return(one.fluxVal < two.fluxVal);}


static void normval_weighted_median(PixVal *Values, const int Npix, 
				    double &Flux, double &FluxWeight)
			       
{
  /* optimize allocation of local array */
  static double* weights = NULL;
  static double* sum_weights = NULL;
  static int weight_size = 0;

  if (Npix == 0)
    {
      Flux = 0;
      FluxWeight = 0;
      return;
    }
     
  if (weight_size < Npix)
    {
      if (weights) {delete[]weights;delete[]sum_weights;}
      weights = new double[Npix];
      sum_weights = new double[Npix];
      weight_size = Npix;
    }
 
  //sort (increasing fluxVal) to compute median
  sort(Values, Values+Npix, PixVal_fluxVal_compare);

  // get a weighted median (i.e sum of weights equal on both sides)
  double sumWeights = 0;
  for (int i=0; i<Npix; ++i)
    {
      const PixVal *value = Values+i;
      /* see comments below for the value of the weight (= 1/Variance) */
      weights[i] = value->totalWeight * sqr(value->photomRatio);
      sum_weights[i] = sumWeights + weights[i];
      sumWeights = sum_weights[i];
    }

  double halfSumWeights = 0.5*sumWeights;
  int k = 0;
  for (; k<Npix; ++k)
    if (sum_weights[k] > halfSumWeights) break;
  
  if (k==0)
    { 
      double v1 = Values[1].fluxVal;
      if (!sumWeights)
	Flux = 0;
      else
	Flux = (Values[0].fluxVal - v1) * weights[0] / sumWeights + v1;
    }
  else
    {
      //weighted average of k-1 and k with resp. weights sum_weight[k-1],[k]
      double v1 = Values[k].fluxVal;
      if (!sumWeights)
	Flux = 0;
      else
	Flux = (Values[k-1].fluxVal - v1) * sum_weights[k-1] / sumWeights + v1;
    }
  FluxWeight = sumWeights/1.5;
 
  /* this factor is not totally full proof: it is approximately valid
     asymptotically. We would like to get a better approximation that 
     holds for this weighted median. 
     Kendal & Stuart propose an approximation for unweighted medians.
  */
}


// does the weighted average of fluxVal
/* 
   by definition : pixVal = fluxVal * photomRatio
   localWeight = 1/Var(pixVal)
   
   hence Var(fluxVal) = 1/(localWeight * photomRatio^2)
   and the weighted (by the inverse of variances) average of fluxVal reads:
   sum(localWeight * photomRatio * pixVal)
   --------------------------------------
   sum(localWeight * photomRatio^2)
   
   we use pixVal rather than fluxVal at the numerator to avoid one division
*/

static void flux_and_weight(const PixVal *Values, const int Npix, 
			    double &Flux, double &WeightFlux)
{
  double sum_flux_num = 0;
  double sum_flux_deno = 0;
  double sum_weight_deno = 0;
  for (int i=0; i<Npix; ++i)
    {
      const PixVal *value = Values+i;
      if (value->keep)
	{
	  double product = value->totalWeight * value->photomRatio; 
	  sum_flux_num += value->pixVal * product;
	  sum_weight_deno += sqr(product) * value->pixVar;
	  sum_flux_deno += product * value->photomRatio;
	}
    } 
  
  if (!sum_flux_deno)
    {
      Flux = 0;
      WeightFlux = 0;
    }
  else
    {
      Flux = sum_flux_num/sum_flux_deno;
      WeightFlux = sqr(sum_flux_deno)/sum_weight_deno;
    }
}

// this routine assumes that 1/localWeight is the variance of rawVal.
static void clipped_weighted_average(PixVal *Values, const int Npix, 
				     const double NSigCut,
				     double &Flux, double &WeightFlux)
{
  double median, median_weight;
  
  normval_weighted_median(Values, Npix, median, median_weight);
  
  // cut around median ...
  double cut2 = sqr(NSigCut);
  for (int i=0; i< Npix; ++i)
    {
      PixVal *value = Values+i;
      // fluxVal = imageVal/photomRatio, and Var(imageVal) = 1/localWeight
      // hence Var(fluxVal) = 1/(localWeight*photomRatio^2)
      if (sqr(value->fluxVal - median) > 
	  cut2 / (value->localWeight * sqr(value->photomRatio)) )
	{
	  value->keep = false;
	}
    }
  // ... and compute average
  flux_and_weight(Values, Npix, Flux, WeightFlux);
}

// continuous thresholding: update the weight function
// weight = weight / ( 1 + weight*resid^2)
// to allow to deweight high residuals
// this routine assumes that 1/localWeight is the variance of rawVal.
static void adaptive_weighted_average(PixVal *Values, const int Npix, 
				      double &Flux, double &WeightFlux)
{
  // initial guess is median
  normval_weighted_median(Values, Npix, Flux, WeightFlux);
  double *origWeight = new double[Npix];
  for (int i=0; i< Npix; ++i) {
    PixVal *value = Values+i;
    origWeight[i] = value->localWeight;
  }
  int iter = 0, maxIter=5;
  do {
    double *pw = origWeight;
    for (int i=0; i< Npix; ++i, ++pw) {
      PixVal *value = Values+i;
      value->localWeight = *pw / (1. + sqr(value->fluxVal-Flux) * value->localWeight);
    }
    flux_and_weight(Values, Npix, Flux, WeightFlux);
  } while (iter++ < maxIter);
  delete [] origWeight;
}

#include "fitsslice.h"
#include <time.h>
 
bool ImageSum::MakeFits()
{

  string fileName = ReducedImage::FitsName();
  if (FileExists(fileName)) return true;
  if (components.size() == 1)
    {
      string name = components.front().Ri->FitsName();
      if (!FileExists(name)) components.front().Ri->MakeFits();
      MakeRelativeLink(components.front().Ri->FitsName().c_str(), FitsName().c_str());
      name = components.front().Ri->FitsWeightName();
      if (!FileExists(name)) components.front().Ri->MakeWeight();
      MakeRelativeLink(components.front().Ri->FitsWeightName().c_str(), FitsWeightName().c_str());
      return true;
    }
  cout <<" ImageSum: entering image stacking for " << Name() << endl;
  clock_t tstart = clock();
  FitsParallelSlices imageSlices(20);
  FitsParallelSlices *weightSlices = NULL;
  // if NoWeightsAtAll, we have to ignore input weights (this is silly, but?)
  // we will use weightSlice to trace wether we use input weightMaps or not
  if (weightingMethod != NoWeightsAtAll)
    weightSlices = new FitsParallelSlices(20);
  
  cout << " ImageSum: stacking using  " 
       << name_of_stackingMethod(stackingMethod) << endl
       << " with weighting scheme : " 
       << name_of_weightingMethod(weightingMethod) << endl;

  PixVal *pixValues = new PixVal[components.size()];
  
  cout << " ImageSum: setting up components" << endl;
  double averageWeightSum = 0;
  for (size_t k=0; k < components.size(); ++k)
    {
      ReducedImage &ri = *(components[k].Ri);
      ri.MakeFits();
      imageSlices.AddFile(ri.FitsName());
      double average_local_weight = 1;
      if (weightSlices) 
	{
	  string weightName = ri.FitsWeightName();
	  //ri.MakeWeight();
	  {
	    /* compute an average image weight for stack attributes 
	       (seeing...) set in FitsHeaderFill */
	    FitsImage weight(weightName);
	    float average,sigma;
	    weight.SkyLevel(&average,&sigma);
	    average_local_weight = average; 
	  }
	  weightSlices->AddFile(weightName);
	}
      else
	{
	  cerr << " ImageSum: do not use weight image for " 
	       << components[k].Ri->Name() << endl ;
	  average_local_weight = 1;
	}
      components[k].averageWeight = average_local_weight
	                            *components[k].globalWeight
                                    *components[k].photomRatio;
      averageWeightSum += components[k].averageWeight;
    }

  // Not sure this is meaningful because there are attributes 
  // which add up, and attributes which average (seeing)
  if (averageWeightSum) for (size_t k=0; k < components.size(); ++k)
      components[k].averageWeight /= averageWeightSum;

  Image stackImage(components[0].Ri->XSize(), components[0].Ri->YSize());
  Image weightImage(stackImage);

  size_t numberOfImages = components.size();
  int nx = stackImage.Nx();
  int ny = stackImage.Ny();

  // set up done : loop on pixels
  cout << " ImageSum: start looping on pixels " << endl;
  do
    {
      for (int j=0; j < imageSlices.SliceSize(); ++j) 
	{

	  int j_image = imageSlices.ImageJ(j);
	  if (j_image >= ny)
	    cout << " ImageSum: CATASTROPHE :j_image >= ny " << j_image 
		 << ' ' << ny << endl; 
	  for (int i=0; i< nx; ++i)
	    { 
	      int npix = 0;
	      for (size_t k = 0; k<numberOfImages; ++k) 
		{
		  /* there are quantities initialized here which do not vary 
		     from pixel to pixel. They are anyway (re) initialised
		     at every pixel because the statistics involved 
		     later operate sorts on the array. 
		     So do not move these initializations from inside the loop
		     to outside the loop */
		  /* on top of that, we ignore here pixels with zero weight */
		  double localWeight=1.;
		  if (weightSlices)
		    localWeight = (*(*weightSlices)[k])(i,j);
		
		  // ignore pixel which have a null local or global weight
		  if (localWeight == 0) continue;
		  if (components[k].globalWeight == 0) continue;

		  // assign "pointer"
		  PixVal &pixStuff = pixValues[npix];

		  // compute pixel Variance
		  if (weightSlices) // localWeight !=0 because of test above
		    pixStuff.pixVar = 1/localWeight; 
		  else pixStuff.pixVar = components[k].backVar;

		  pixStuff.globalWeight = components[k].globalWeight;
		  pixStuff.pixVal = (*imageSlices[k])(i,j);
		  pixStuff.localWeight = localWeight;
		  pixStuff.totalWeight = localWeight * pixStuff.globalWeight;
		  pixStuff.photomRatio = components[k].photomRatio;
		  pixStuff.fluxVal = pixStuff.pixVal/pixStuff.photomRatio;
		  pixStuff.keep = true;
		  npix++;
		}
	      // data collected : carry out the statistics
	      double stackVal; // resulting flux estimate
	      double weightVal = 0; // its weight
	      if (stackingMethod == ClippedWeightedAverage)
		{
		  clipped_weighted_average(pixValues, npix, 5.,
					   stackVal, weightVal);
		}		  
	      else if (stackingMethod == WeightedAverage)
		{
		  flux_and_weight(pixValues, npix, stackVal, weightVal);
		}
	      else if (stackingMethod == Median)
		{
		  normval_weighted_median(pixValues, npix, 
					  stackVal, weightVal);
		}
	      else if (stackingMethod == AdaptiveWeightedAverage)
		{
		  adaptive_weighted_average(pixValues, npix, stackVal, weightVal);
		}
	      stackImage(i, j_image) = stackVal;
	      weightImage(i, j_image) = weightVal;
	    }
	}
      
      if (weightSlices)  weightSlices->LoadNextSlice();
    } while (imageSlices.LoadNextSlice());
  
  if (weightSlices) delete weightSlices;
  delete [] pixValues;

  clock_t tend = clock();
  cout << " ImageSum: CPU for stacking " 
       << float(tend-tstart)/float(CLOCKS_PER_SEC) << endl;
  {
    FitsHeader headRef(components[0].Ri->FitsName());
    FitsImage stackFits(fileName, headRef, stackImage);
    stackFits.SetWriteAsFloat();
    FitsImage weightFits(FitsWeightName(), headRef, weightImage);
    // to make sure that 0 remain 0 after FITS write/read in 16 bits format:
    weightFits.PreserveZeros(); 
  }
  FitsHeaderFill();
  cout << " ImageSum: " << Name() << " Image/Weight stat : " 
       << ImageAndWeightError(stackImage, weightImage) << endl;
  return true;
}

bool ImageSum::MakeWeight()
{
  // weight image is constructed in parallel to image.
  if (!FileExists(FitsWeightName())) return MakeFits();
  else return true;
}

void ImageSum::FitsHeaderFill()
{
  double seeing = 0;
  double backLevel = 0;
  double backVar = 0;
  double saturation = 0;
  double originalsatur = 0;
  double exposure = 0;
  double readnoise = 0;
  double flatnoise = 0;
  double julianDate = 0;
  double airmass = 0;
  double originalskylevel = 0;
  double num_gain = 0;
  double deno_gain = 0;
  // string componentsNames;

  for (size_t k=0; k<components.size(); ++k)
    {
      cout << " ImageSum: Component #" << k << endl;
      Component &component = components[k];
      component.dump();
      ReducedImage* redIm = component.Ri;

      //componentsNames += (redIm->Name() + " ");

      /* SEEING */
      seeing += component.averageWeight*redIm->Seeing();
  
      /* BACKGROUND */
      // by construction, it should be zero for both input and output.
      backLevel += component.averageWeight*redIm->BackLevel();

      /* BACKGROUND Variance */
      backVar += sqr(component.averageWeight)*component.backVar;
  

      /* SATURATION */
      saturation += component.averageWeight * redIm->Saturation();
  
      /* EXPOSURE */
      exposure += redIm->Exposure();

      /* DATE */
      julianDate += component.averageWeight * redIm->JulianDate();
  
      /* AIRMASS */
      airmass += component.averageWeight * redIm->Airmass();

      /* READ OUT NOISE*/
      readnoise += sqr(component.averageWeight * redIm->ReadoutNoise());
  
      /* FLAT NOISE */
      flatnoise += sqr (component.averageWeight * redIm->FlatFieldNoise()
			*redIm->OriginalSkyLevel());
  
      /* ORIGINAL SKY LEVEL */
      originalskylevel += component.averageWeight * redIm->OriginalSkyLevel();

      /* ORIGINAL BACKGROUND */
      originalsatur += component.averageWeight * redIm->OriginalSaturation();
  

      /* gain */
      double gain = redIm->Gain();
      if (gain!= 0) 
	/* I am afraid this formula only holds for the OptimalWeighting,
	   WeightedAverage case. */
	{
	  deno_gain += sqr(component.averageWeight)*component.photomRatio/gain;
	  num_gain += component.averageWeight*component.photomRatio;
	}
    }

  // Upadte the usable part of the image
  //  SetUsablePart(intersection);
  cout << " ImageSum: setting the header " << endl;
  SetSeeing(seeing);
  SetBackLevel(backLevel);
  double sigmaBack = sqrt(backVar);
  SetSigmaBack(sigmaBack);
  float average, sigma;
  {
    FitsImage image(FitsName());
    image.SkyLevel(&average,&sigma);
  }
  SetNoisePow(sigma);
  SetJulianDate(julianDate);
  SetAirmass(airmass);
  SetSaturation(saturation);
  SetExposure(exposure);
  readnoise = sqrt(readnoise);
  SetReadoutNoise(readnoise);

  SetOriginalSkyLevel(originalskylevel);
  flatnoise = sqrt(flatnoise)/originalskylevel;
  SetFlatFieldNoise(flatnoise);

  SetOriginalSaturation(originalsatur);

  if (deno_gain != 0) SetGain(sqr(num_gain)/deno_gain);
  SetPhotomReference(photomReferenceName);

  FitsHeader head(FitsName(), RW);
  head.AddOrModKey("STACKMET", name_of_stackingMethod(stackingMethod), 
		   " ImageSum param ");
  head.AddOrModKey("WEIGHMET", name_of_weightingMethod(weightingMethod), 
		   " ImageSum param ");

  // final frame can be large
  if (head.HasKey("PKAUNION")) {
    head.RmKey("XNEWBEG");
    head.RmKey("YNEWBEG");
    head.RmKey("XNEWEND");
    head.RmKey("YNEWEND");
  }

   if (zero_point_ref>0)
    {
      // par construction, sum alignee sur ref photom
      cout << " ImageSum: zero point for stack " << Name() << " : " 
	   << zero_point_ref << endl ;
      SetZZZeroP(zero_point_ref, "computed from ref for stack");
    }


  //  head.AddOrModKey("STACKCMP", componentsNames, 
  //		   " DbImage names of components");

  // add history lines instead of a key. It crashed with long names... 
  // I suspect the string capacity holding with g++<3.2
  head.AddHistoryLine("DbImage names of components");
  for (size_t k=0; k<components.size(); ++k)
    head.AddHistoryLine(components[k].Ri->Name());
    
}

ReducedImageList ImageSum::Components() const
{
  ReducedImageList list;
  for (ComponentCIterator i = components.begin(); i != components.end(); ++i)
    list.push_back(i->Ri);
  return list;
}

bool ImageSum::MakeCatalog() 
{
  if (FileExists(CatalogName()))  return true;
  if (components.size() == 1)
    {
      MakeRelativeLink(components.front().Ri->CatalogName().c_str(),
		       CatalogName().c_str());
      return true;
    }
  // mode overwrite, ne fait pas la carte des etoiles saturees,
  // ne soustrait pas le fond, utilise le sigma du header.
  //bool ok = ReducedImage::MakeCatalog_ImageBizarre();
  bool ok = MakeCatalog_ImageBizarre(); 
  //appellera MakeWeight de imagesum car virtual
  
  if( !HasSatur())
    MakeSatur();
  if (HasSatur())// flagger etoiles saturees
    {
      FitsImage satur(FitsSaturName());
      SEStarList stl(CatalogName());
      FlagSatFromImage(stl, satur);
      stl.write(CatalogName());
    }
  return(ok); 
}  

bool ImageSum::MakeAperCat()
{
  if (FileExists(AperCatalogName())) return true;
  if (components.size() == 1)
    {      
      string name = components.front().Ri->AperCatalogName();
      if (!FileExists(name)) components.front().Ri->MakeAperCat();
      MakeRelativeLink(components.front().Ri->AperCatalogName().c_str(),
		       AperCatalogName().c_str());
      return true;
    }  
  return ReducedImage::MakeAperCat();
}  

bool ImageSum::MakeStarCat()
{
  if (FileExists(StarCatalogName())) return true;
  if (components.size() == 1)
    {      
      string name = components.front().Ri->StarCatalogName();
      if (!FileExists(name)) components.front().Ri->MakeStarCat();
      MakeRelativeLink(components.front().Ri->StarCatalogName().c_str(),
		       StarCatalogName().c_str());
      return true;
    }
  return ReducedImage::MakeStarCat();
}  

bool ImageSum::MakeDead()
{
  if (FileExists(FitsDeadName())) return true;
  cout << " ImageSun: making Dead image for " << Name() << endl;
  if (components.size() == 1)
    {
      MakeRelativeLink(components.front().Ri->FitsDeadName().c_str(),
		       FitsDeadName().c_str());
      return true;
    }
  ReducedImageList components = Components();
  return BoolImageOr(components, &ReducedImage::FitsDeadName, 
		     &ReducedImage::MakeDead, FitsDeadName());
}

bool ImageSum::MakeSatur()
{
  if (FileExists(FitsSaturName())) return true;
  cout << " ImageSum: making Satur image for " << Name() << endl;
  if (components.size() == 1)
    {
      string name = components.front().Ri->FitsSaturName();
      if ( !FileExists(name) )
	components.front().Ri->MakeSatur();
      MakeRelativeLink(components.front().Ri->FitsSaturName().c_str(),
		       FitsSaturName().c_str());
      return true;
    }
  ReducedImageList components = Components();
  for (ReducedImageIterator i=components.begin(); i != components.end(); ++i)  
    {
      string name = (*i)->FitsSaturName();
      if (!FileExists(name) )
	(*i)->MakeSatur();
    }
  return BoolImageOr(components, &ReducedImage::FitsSaturName, 
		     &ReducedImage::MakeSatur, FitsSaturName());
}


/*****************************************************************************/
//! Align (on Reference) and sum images. ToDo can be constructed using DoFits 
//DoCatalog DoDead DoSatur. 
ImageSum* ImagesAlignAndSum(const ReducedImageList &ToSum, 
			    const ReducedImage &Reference, 
			    const string &SumName, const int ToDo,
			    const ReducedImage *PhotomReference,
			    const WeightingMethod AWMethod, 
			    const StackingMethod ASMethod,
			    const PhotoRatioMethod APMethod)

{
  ReducedImageList transformedImages;
  // align
  int count = ImagesAlign(ToSum, Reference, transformedImages, ToDo);
  //sum
  if (count == 0) return NULL;
  ImageSum *sum;
  sum = new ImageSum(SumName, transformedImages, PhotomReference, AWMethod, ASMethod, APMethod);
  sum->Execute(DoFits | ToDo);
  
  return sum;
}



bool ImagesAlignAndSum(const vector<string> &ToSum, 
		       const string &Reference, 
		       const string &SumName, const int ToDo)
{
  ReducedImage geomRef(Reference);
  ReducedImageList reducedImages;
  if (!geomRef.IsValid())
    {
      cerr << " ImagesAlignAndSum : was " << Reference 
	   << " actually produced ?? " << endl;
      return false;
    }
  for (size_t i=0; i<ToSum.size() ; ++i)
    {
      /* I guess I should refer to the ReducedImage via the virtual 
	 constructor (ReducedImageNew) rather than the real one. 
	 The problem is that the virtual one does not actually work as long as
	 all the data inside the ReducedImage is not written to disk 
	 which is not done. */
      string currentName = ToSum[i];
      ReducedImage *current = new ReducedImage(currentName);
      if (!current->IsValid())
	{
	  cerr << " ImagesAlignAndSum : was " << currentName 
	       << " actually produced ?? " << endl;
	  continue;
        }
      reducedImages.push_back(current);
    }
  ReducedImage *result = ImagesAlignAndSum(reducedImages, Reference, 
					   SumName, ToDo);
  if (result) delete result;
  else {
    cerr << " ImagesAlignAndSum: failure" << endl;
    return false; 
  }
  return true;
}



