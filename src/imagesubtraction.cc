#include <algorithm>

#include "imagesubtraction.h"
#include "kernelfit.h"
#include "fileutils.h"
#include "fitsimage.h"
#include "imageback.h"
#include "gtransfo.h"


#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif

string SubtractedName(const string &RefName, const string &NewName) {return NewName+'-'+RefName;}

ImageSubtraction::ImageSubtraction(const string &Name,  
				   const ReducedImageRef RefImage,  const ReducedImageRef NewImage,
				   const bool NoSwap ) : 
  ReducedImage(Name), KernelFitter(RefImage, NewImage, NoSwap) , Ref(RefImage), New(NewImage)
{
  if (!FileExists(Name))
    Create("here");
}

ImageSubtraction::ImageSubtraction(const string &Name,  const ReducedImageRef  RefImage,  const ReducedImageRef NewImage, 
				   const KernelFitter &APreviousFit) :
  ReducedImage(Name), KernelFitter(APreviousFit), Ref(RefImage), New(NewImage)
{
  Create("here");
}

ReducedImage *ImageSubtraction::Clone() const
{
  /* this routine should be useless. I don't know why it is here ... */
  /* 
  ReducedImage *clone = new ImageSubtraction(*this);
  return clone;
  */
  return NULL;
}


#ifdef STORAGE
// Code to generate Delphine's plots (subtraction quality checking) , as they are called.

#include "quali_box.h"

// this routine creates a kind of ntuple that allows to study 
// the subtraction quality  (the famous Delphine's plot)
static void  quali_plots(const ReducedImage *Ref, const ReducedImage *New, 
	    const FitsImage *subImage, const string LogName)
{
  const Image *imArray[3];
  FitsImage refImage(Ref->FitsName());
  FitsImage newImage(New->FitsName());
  imArray[0] = &refImage;
  imArray[1] = &newImage;
  imArray[2] = subImage;
  double fond[3] = {Ref->BackLevel(), New->BackLevel(), 0};  /* the last one is zero because the kernelfit is supposed to achieve that ... */
  double sigmaFond[3] = {Ref->SigmaBack(), New->SigmaBack(),
			 sqrt(sqr(Ref->SigmaBack())+sqr(New->SigmaBack()))} ;
  /* bon d'accord, le dernier est ``faux'' mais chez Delphine aussi. 
     alors comme c'est pour faire des comparaisons... */

  SEStarList refList(Ref->CatalogName());
  SEStarList newList(New->CatalogName());
  ofstream qual_log(LogName.c_str());
  GtransfoIdentity ident;

  Test_Qual_N_same(3, *subImage, refList, newList,
		   &ident, imArray, fond, sigmaFond, 10, qual_log);
  qual_log.close();
}                 
#endif /* STORAGE */


static double sqr(double x) { return x*x;}


bool ImageSubtraction::MakeFits()
{
  string fileName = FitsName();
  if (FileExists(fileName)) return true;
  if (!Ref->MakeFits() ||  !New->MakeFits())
    {
      cerr << " ERROR : for ImageSubtraction " << Name() << endl
	   << " Could not make images of (both) subtraction terms :" << endl
	   << Ref->Name() << " and " << New->Name() << endl;
      return false;
    }
  // fit the convolution kernel
  DoTheFit();
  // build the subtracted image
  FitsHeader href(Ref->FitsName());
  FitsImage subImage(fileName, href);
  Image &theSubtraction = subImage;
  /* by convention the subtraction photometric scale is the same as the ref.
     - photomRatio is New/Ref
     - Important point : the same trick has to be applied to the variance computation. 
     use subblocks in the following code to get rid ASAP of temporary images
  */
  cout << " Subtracting " << New->Name() << " - " << Ref->Name() << endl;
  double photomRatio = KernAtCenterSum();
  if (RefIsBest())
    {
      {
	FitsImage worst(New->FitsName());
	theSubtraction = worst;
      }
      {
	FitsImage best(Ref->FitsName());
	Image convBest(best.Nx(), best.Ny());
	ImageConvolve(best,convBest);
	theSubtraction -= convBest;
      }
      double factor = 1./photomRatio; // "*" is far faster than "/"
      theSubtraction *=factor;
    }
  else 
    {
      {
	FitsImage best(New->FitsName());
	ImageConvolve(best,theSubtraction);
      }
      {
	FitsImage worst(New->FitsName());
	theSubtraction -= worst;
      }
    }
  // force floating point values (not sure it is useful)
  subImage.SetWriteAsFloat();
  // update Frame limits
  CommonFrame().WriteInHeader(subImage);
  // link PSF if any
  string worstpsfname = Worst()->ImagePsfName();
  if (FileExists(worstpsfname)) MakeRelativeLink(worstpsfname.c_str(), ImagePsfName().c_str());
  // set some score values
  SetSeeing(LargestSeeing());
  SetBackLevel(0.0); // by construction
  SetSaturation(Worst()->Saturation()," worst image saturation");
  double sigmaBack = sqrt(sqr(Best()->SigmaBack()*photomRatio) + sqr(Worst()->SigmaBack()));
  SetSigmaBack(sigmaBack);
  double readnoise = sqrt(sqr(Worst()->ReadoutNoise())+sqr(photomRatio*Best()->ReadoutNoise()));
  SetReadoutNoise(readnoise);
  /* originaly , we added some info on the origin of this image. Not sure it is useful */
  /*
  AddOrModKey("KERNREF",Best()->Name().c_str(), " name of the seeing reference image");
  AddOrModKey("KERNCHI2",Chi2(), " chi2/dof of the kernel fit");
  AddOrModKey("PHORATIO",photomRatio, " photometric ratio with KERNREF"); 
  */

  //  quali_plots(Ref, New, &subImage, Dir()+"/qual.log");

  double zero_point = Ref->AnyZeroPoint();
  cerr << "Zero Point as taken from reference stack:" 
       << zero_point << " propagated to sub. " << endl ;
  SetZZZeroP(zero_point, "as taken from reference stack");  

  // subtract an ImageBack
  FitsImage *pweight = NULL;
  
  if (MakeWeight())
    {
      pweight = new FitsImage(FitsWeightName());
      cout << " using weights : " << FitsWeightName() << endl;
    }
  int backMesh = 32;
  ImageBack b(theSubtraction, backMesh, pweight);
  if (pweight) delete pweight;

  Image *back = b.BackgroundImage();
  theSubtraction -= *back ;
  delete back;
  subImage.AddOrModKey("BACK_SUB", true, 
		       pweight? " subtracted weighed background" :
		      " subtracted unweighted background");
  subImage.AddOrModKey("BACKMESH", backMesh, 
		       " mesh size used for back computation " );
  return true;
}

#include "convolution.h"



  /* this routine computes the weight map of the subtraction:
     - weights from best image
     - account for satur in this image
     - go to variances
     - convolve (by the squared kernel form the image kernel fit)
     - take weights from worst and account for satur in worst.
     - combine both maps to get total variance
     - normalize with the same convention as the image
     - go back to weights
     - zero the sides of the weight map
     - done. easy isn't it?
  */

bool ImageSubtraction::MakeWeight()
{
  string fileName = FitsWeightName();
  if (FileExists(fileName)) return true;
  if (!MakeFits()) return false; // Hard way to get the convolution kernel.

  cout << " making WeightImage  " << fileName << endl;
  FitsHeader imageHeader(FitsName()); 

  FitsImage  weights(FitsWeightName(), imageHeader);
  Image &weightImage = weights; // handle to resolve FitsHeader/Image ambiguities
  { // to open temporarily Best weight image.
    FitsImage best_weight(Best()->FitsWeightName());
    weightImage = best_weight;
  }
  
  // account for satur pixels in Best
  if (Best()->MakeSatur())
    {
      FitsImage satur(Best()->FitsSaturName());
      weightImage *= (1 - satur);
    }
  else
    {
      cerr << " When making weight map for "<< Name() 
	   << ", cannot get satur from " << Best()->Name() << endl;
    }

  /* add a small constant to weights so that variances of 
     zero weight pixels remain finite */
  Pixel minVal,maxVal;
  weights.MinMaxValue(&minVal, &maxVal);
  double eps = maxVal*1e-10;
  weights += eps;

  // now go to variances
  Pixel *pend = weights.end();
  for (Pixel *p = weights.begin(); p < pend; ++p) (*p) = 1./(*p);

  /* convolve with squared image kernel */
  cout << " convolving variance " << endl;
  {
    Image weightsCopy(weights);
    VarianceConvolve(weightsCopy, weights);
  }

  // load the other weight map
  FitsImage weights_worst(Worst()->FitsWeightName());

  // zero saturated pixels in Worst weight map
  if (Worst()->MakeSatur())
    {
      FitsImage satur(Worst()->FitsSaturName());
      weights_worst *= (1 - satur);
    }
  else
    {
      cerr << " When making weight map for "<< Name() 
	   << ", cannot get satur from " << Worst()->Name() << endl;
    }

  Pixel threshold = eps*100; // value under which we go to zero.
  cout << " threshold under which weight = 0 " << threshold << endl;
  pend = weights.end();
  Pixel *pworst = weights_worst.begin();
  for (Pixel *p = weights.begin(); p < pend; ++p, ++pworst)
    { // the following expression works for *pworst == 0
      (*p) = (*pworst)/((*p)*(*pworst)+1.);
      if (*p < threshold) { (*p) = 0;}
    }

  // weightImage is now an actual weight.

  /* account for the fact that the subtraction is by convention
     expressed in units of the ref. The weight map we have is
     correctly normalized only if New was convolved (and
     photometrically matched to ref by the kernel fit).  Apply the same
     factor as the one applied to the image (however squared and
     inverted) in PsfMatch::Subtraction() */

  if (RefIsBest())
    {
      weightImage *= sqr(KernAtCenterSum());
    }


  /* set to 0 the "side bands" where the convolution did not go */


  Frame aframe(weightImage); // whole image

  // size of the dead band due to variance convolution 
  Kernel kern_at_center; // get a kernel to grab its size.
  KernAllocateAndCompute(kern_at_center, weights.Nx()/2., weights.Ny()/2);

  int bandx = kern_at_center.Nx()/2 -1;
  int bandy = kern_at_center.Ny()/2 -1;
  // shrink the frame :
  aframe.CutMargin(bandx, bandy); 
  // outside of this, weights are to be set to 0
  
  weights.Masking(aframe, 0.);
  cout << " masking frame for weights " << aframe << endl;

  weights.PreserveZeros(); // 0 remain 0 on R/W operations
  weights.ModKey("BITPIX",16); // 16 bits are enough
  return true;
}


bool ImageSubtraction::MakeDead()
{
  if (FileExists(FitsDeadName())) return true;
  ReducedImageList list(false);
  list.push_back(Ref);
  list.push_back(New);
  bool return_value = BoolImageOr(list, 
				  &ReducedImage::FitsDeadName, 
				  &ReducedImage::MakeDead, FitsDeadName());
  FitsImage dead(FitsDeadName(), RW);
  cout << " the frame for dead of " << Name() << " is " << UsablePart() << endl;
  cout << " # of dead pixels before " << dead.SumPixels() << endl;
  dead.Masking(UsablePart(),1);
  cout << " # of dead pixels after " << dead.SumPixels() << endl;
  return return_value;
}

bool ImageSubtraction::MakeSatur()
{
  if (FileExists(FitsSaturName())) return true;
  ReducedImageList list(false);
  list.push_back(Ref);
  list.push_back(New);
  return BoolImageOr(list, &ReducedImage::FitsSaturName, &ReducedImage::MakeSatur, FitsSaturName());
}

bool ImageSubtraction::MaskSatur()
{
   if (!FileExists(FitsName()) || !FileExists(FitsSaturName())) return false; 
   FitsImage img(FitsName(),RW);
   FitsImage satur(FitsSaturName());
   Image &i = img;
   i *= (1.-satur);
   return true;
 }

bool ImageSubtraction::MaskNullWeight()
{
 if (!FileExists(FitsName()) || !FileExists(FitsWeightName())) return false;
 FitsImage img(FitsName(),RW);
 FitsImage weight(FitsWeightName());
 weight.Simplify(1e-30);
 Image &i = img;
 i *= (weight);
 return true;
}



string ImageSubtraction::CandName() const
{
  return Dir()+"/cand.list" ;
}

string ImageSubtraction::AllCandName() const
{
  return Dir()+"/allcand.list" ;
}

string ImageSubtraction::CandCutName() const
{
  return Dir()+"/cand_cut.list" ;
}

string ImageSubtraction::CandScanName() const
{
  return Dir()+"/cand_scan.list" ;
}

string ImageSubtraction::CandCutScanName() const
{
  return Dir()+"/cand_cut_scan.list" ;
}


#include "toadscards.h"



#include "detection.h"

bool ImageSubtraction::RunDetection(DetectionList &Detections,
				    const BaseStarList* Positions, string name,bool fix_pos)
{
  if (name == "") name =DetectionsName();
  if (!MakeFits() || !MakeWeight())
    {
      cerr << " cannot make (sub) catalog for " << Name() 
	   << " without both image and weights " << endl;
      return false;
    }
  if (!Positions && FileExists(DetectionsName()))
    {
      Detections.read(DetectionsName());
      return true;
    }

  // compute some informative statistics 
  {
    FitsImage im(FitsName());
    FitsImage w(FitsWeightName());
    cout << Name() << " Image/Weight stat : " << ImageAndWeightError(im,w) << endl;
  }

  // filter size is equal to seeing until we refine it

  if (!Positions)
    {
      {//block to save memory
	DetectionProcess detectionProcess(FitsName(), FitsWeightName(), 
					  Seeing(), Seeing());
	detectionProcess.DoDetection(Detections);
      }
      {//block to save memory
	/* fill in the Detection's block that concerns the ref:
	   - flux on the ref under the SN (measured in the same conditions
	   as the SN (this is why we use the sub Seeing() 
	   rather the ref.Seeing())
	   - nearest object
	*/
	DetectionProcess refDet(Ref->FitsName(), Ref->FitsWeightName(), 
				Seeing(), Seeing());
	refDet.SetScoresFromRef(Detections, *Ref);
      }
      Detections.write(name);
    }
  else
    { 
      {
      DetectionProcess detectionProcess(FitsName(), FitsWeightName(), 
					  Seeing(), Seeing());
      detectionProcess.DetectionScoresFromPositions(*Positions, Detections,fix_pos);
      }
      {//block to save memory
	/* fill in the Detection's block that concerns the ref:
	   - flux on the ref under the SN (measured in the same conditions
	   as the SN (this is why we use the sub Seeing() 
	   rather tha ref.Seeing())
	   - nearest object
	*/
	DetectionProcess refDet(Ref->FitsName(), Ref->FitsWeightName(), 
				Seeing(), Seeing());
	refDet.SetScoresFromRef(Detections, *Ref);
       
      }
      Detections.write(name);
    }
  cout << "zero the pixel with null weight" << endl;
  MaskNullWeight();
  
  return true;
}


bool ImageSubtraction::MakeCatalog()
{
  if (FileExists(DetectionsName())) return true;
  DetectionList detections;
  return RunDetection(detections);
}


#ifdef USE_ROOT
/* 
RUN_ROOTCINT
LINKDEF_CONTENT : #pragma link C++ class PsfMatch+;
LINKDEF_CONTENT : #pragma link C++ class ImageSubtraction+;

*/
ClassImp(ImageSubtraction);

#include "root_dict/imagesubtractiondict.cc"

#endif /* USE_ROOT */
