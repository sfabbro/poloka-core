#include <fstream>

#include "reducedutils.h"
#include "starmatch.h"
#include "listmatch.h"
#include "apersestar.h"
#include "imageutils.h"
#include "vutils.h" // for DArrayMedian
#include "transformedimage.h"
#include "fitsimage.h"
#include "wcsutils.h"
#include "toadscards.h"
#include "allreducedimage.h"
#include "polokaexception.h"

void LoadForMatch(const ReducedImage& Im, BaseStarList& BList, const double& MinSN) {
  
  SEStarList sList;
  if (Im.HasAperCatalog()) {
    AperSEStarList asList(Im.AperCatalogName());
    sList = AperSE2SE(asList);
  } else if (Im.HasCatalog()) {
    sList.read(Im.CatalogName());
  }

  for (SEStarIterator it=sList.begin(); it != sList.end(); ) {
    SEStar *pstar = *it ;
    // not in bad zone, keep satur, only with better 0.3 pixels accuracy
    if (pstar->FlagBad() == 0 && pstar->Flag() < 8 && 
	pstar->flux > MinSN * pstar->EFlux())
      //pstar->vx + pstar->vy < 0.1)
      ++it;
    else 
      it = sList.erase(it);
  }  
  BList = *SE2Base(&sList);
}

static bool ListMatchCheck(const BaseStarList& Src, const BaseStarList& Dest,
			   const GtransfoRef Transfo, const double& DistMax, const double& MinMatch=0.1) {
    
  if (!Transfo) return false;

  StarMatchList *match = ListMatchCollect(Src, Dest, Transfo, DistMax);
  if (!match) return false;
  if (match->empty()) {
    cerr << " ListMatchCheck: empty match\n";
    delete match;
    return false;
  }
  size_t nmin = 100;
  nmin = min(size_t(min(Src.size(), Dest.size()) * MinMatch), nmin);
  cout << " ListMatchCheck: min stars to match = " << nmin << endl;
  size_t ngood = StarMatchCheck(match);
  delete match;
  return (ngood > nmin);
}

GtransfoRef FindTransfoFromWCS(const ReducedImage& Src, const ReducedImage& Dest) {

  GtransfoRef src2dest;
  FitsHeader hsrc(Src.FitsName());
  if (!HasLinWCS(hsrc)) return src2dest;
  FitsHeader hdest(Dest.FitsName());
  if (!HasLinWCS(hdest)) return src2dest;
  GtransfoRef wsrc = WCSFromHeader(hsrc);
  GtransfoRef wdest = WCSFromHeader(hdest);
  
  if (wsrc && wdest)
    src2dest = GtransfoCompose(wdest->InverseTransfo(0.001, Dest.UsablePart()), wsrc);

  return src2dest;
}


GtransfoRef FindTransfoCombinatorial(const ReducedImage& Src, const ReducedImage& Dest) {

  BaseStarList srcList;  LoadForMatch(Src, srcList);
  BaseStarList destList; LoadForMatch(Dest, destList);

  MatchConditions cond;
  cond.SizeRatio = Src.PixelSize() / Dest.PixelSize();
  cond.DeltaSizeRatio = 0.1 * cond.SizeRatio;

  // hack for big frames
  if (Src.XSize() > 10000 && Src.YSize() > 10000) {
    cond.NStarsL1 = 2000;
    cond.NStarsL2 = 2000;
  }
  // change conditions with file
  cond.read(DefaultDatacards());
  cond.MaxDist = 4;

  cout << " FindTransfoCominatorial: trying with a combinatorial match\n";
  GtransfoRef transfo = ListMatchCombinatorial(srcList, destList, cond);
  if (transfo) {
    cond.MaxDist = 2;
    cout << " FindTransfoCombinatorial: refining transfo\n";
    transfo = ListMatchRefine(srcList, destList, transfo, cond);
  } else {
    cout << " FindTransfoCombinatorial: no transfo found\n";
    transfo = GtransfoRef();
  }

  if (dynamic_cast<GtransfoIdentity*>((Gtransfo*)transfo) || dynamic_cast<GtransfoPoly*>((Gtransfo*)transfo))
    transfo->Write(GtransfoName(Dest, Src));

  return transfo;
}


GtransfoRef FindTransfo(const BaseStarList& SrcList, const BaseStarList& DestList,
			const ReducedImage& Src, const ReducedImage& Dest) {
  
  cout << " FindTransfo: " 
       << Src.Name() << " (" << SrcList.size() << " stars) to " 
       << Dest.Name() << " (" << DestList.size() << " stars)\n";

  MatchConditions cond;
  cond.SizeRatio = Src.PixelSize() / Dest.PixelSize();
  cond.DeltaSizeRatio = 0.1 * cond.SizeRatio;

  // hack for big frames
  if (Src.XSize() > 10000 && Src.YSize() > 10000) {
    cond.NStarsL1 = 2000;
    cond.NStarsL2 = 2000;
  }
  // change conditions with file
  cond.read(DefaultDatacards());
  
  GtransfoRef transfo;
  
  cout << " FindTransfo: trying with WCS's\n";
  GtransfoRef transfoWcs = FindTransfoFromWCS(Src, Dest);
  if (transfoWcs && !ListMatchCheck(SrcList, DestList, transfoWcs, cond.MaxDist, cond.MinMatchRatio))
    cout << " FindTransfo: WCS not good enough\n";
  else
    transfo = transfoWcs;

  if (!transfo) {
    cond.MaxDist = 4;
    cout << " FindTransfo: trying with a combinatorial match\n";
    transfo = ListMatchCombinatorial(SrcList, DestList, cond);
  }

  if (transfo) {
    cond.MaxDist = 2;
    cout << " FindTransfo: refining transfo\n";
    transfo = ListMatchRefine(SrcList, DestList, transfo, cond);
  } else {
    cout << " FindTransfo: no transfo found\n";
    transfo = GtransfoRef();
  }

  return transfo;
}

GtransfoRef FindTransfo(const ReducedImage& Src, const ReducedImage& Dest) {

  cout << " FindTransfo: " << Src.Name() << " to " << Dest.Name();

  if (Src == Dest) {
    cout << " found identity\n";
    return new GtransfoIdentity();
  }

  string transfoFileName = GtransfoName(Dest, Src);

  if (FileExists(transfoFileName)) {
    cout << " found an existing transfo\n";
    return GtransfoRead(transfoFileName);
  }

  GtransfoRef transfoWcs = FindTransfoFromWCS(Src, Dest);
  if (FileExists(transfoFileName+".wcsonly")) {
    cout << " will use match between WCS\n";
    return transfoWcs;
  }

  if (transfoWcs)
    cout << " found a WCS match\n";
  else
    cout << endl;

  MatchConditions cond;
  cond.SizeRatio = Src.PixelSize() / Dest.PixelSize();
  cond.DeltaSizeRatio = 0.1 * cond.SizeRatio;

  // hack for big frames
  if (Src.XSize() > 10000 && Src.YSize() > 10000) {
    cond.NStarsL1 = 2000;
    cond.NStarsL2 = 2000;
  }
  // change match conditions with user values last
  cond.read(DefaultDatacards());

  BaseStarList srcList;  LoadForMatch(Src, srcList);
  BaseStarList destList; LoadForMatch(Dest, destList);

  GtransfoRef transfo;

  if (transfoWcs && ListMatchCheck(srcList, destList, transfoWcs, cond.MaxDist, cond.MinMatchRatio)) {
    cout << " FindTransfo: refining WCS composition\n";
    transfo = ListMatchRefine(srcList, destList, transfoWcs, cond);
  } else {
    cout << " FindTransfo: search for combinatorial match\n";
    transfo = ListMatch(srcList, destList, cond);
  }

  // only write existing writing routines
  if (dynamic_cast<GtransfoIdentity*>((Gtransfo*)transfo) || dynamic_cast<GtransfoPoly*>((Gtransfo*)transfo))
    transfo->Write(transfoFileName);

  // just touch the file if wcs is the best one
  if (transfo == transfoWcs)
    ofstream t((transfoFileName+".wcsonly").c_str());

  return transfo;
}

string ImageSetName(const ReducedImage& AnImage)
{
  //re-string the date
  int dd,mm,yy;
  sscanf(AnImage.Date().c_str(),"%d/%d/%d",&dd,&mm,&yy);
  char datechr[11];
  sprintf(datechr,"%04d%02d%02d",yy,mm,dd);
  return AnImage.Telescope() + "_" + string(datechr)+ "_" + AnImage.Band();
}

void ArrangeByInstBandNight(const ReducedImageList &Images, vector<ReducedImageList> &ImageSets)
{

  cout << " ArrangeByInstBandNight() : " << Images.size() << " images to arrange " << endl;

  ReducedImageList imlist(false);
  for (ReducedImageCIterator ri = Images.begin() ; ri!=Images.end(); ++ri) 
    imlist.push_back(*ri);

  // image by image empty the imageslist
  while (!imlist.empty())
    {
      // get the first one to make a new NightSet
      ReducedImage *first = imlist.front();
      ReducedImageList thisSet;
      thisSet.push_back(first);

      imlist.erase(imlist.begin());

      // check all others, copy and add entry if from same set delete otherwise
      for (ReducedImageIterator ri = imlist.begin() ; ri!=imlist.end(); )
	{
	  ReducedImage *current = *ri;
	  if (current->SameChipFilterInstNight(*first))
	    {
	      thisSet.push_back(current->Clone()); 
	      ri = imlist.erase(ri);
	    }
	  else ++ri;
	}
      ImageSets.push_back(thisSet);
    }
   cout << " ArrangeByInstBandNight() : arranged " << ImageSets.size() << " image sets" << endl;
}

// to find a good geometric reference
const ReducedImage* BestResolutionReference(const ReducedImageList &ImList)
{
  const ReducedImage *geomRef = ImList.front();
  for (ReducedImageCIterator it=ImList.begin(); it != ImList.end(); ++it)
    {      
      const ReducedImage *current = *it;
      double pixratio = current->PixelSize() / geomRef->PixelSize();
      // if pixel size is 10% smaller, keep it. Otherwise, take best seeing
      if (pixratio < 0.9) geomRef = current;
      else if (current->Seeing() <= geomRef->Seeing() && pixratio < 1.01) 
	geomRef = current;
    }
  cout << " BestResolutionReference() : Found " 
       << geomRef->Name() << ": seeing " << geomRef->Seeing()
       << "  pixel " << geomRef->PixelSize() << endl;
  
  return geomRef;
}

// to get rid of images non overlapping in a list of images
void FilterWithOverlap(const ReducedImageList &Images, const ReducedImage &aReference, 
		       ReducedImageList &overlapImages, const double MinOverlap)
{

  cout << " FilterWithOverlap() : " << Images.size() << " images to check with a min. overlap of " 
       << MinOverlap << " arcmin2" << endl;

  for (ReducedImageCIterator ri = Images.begin(); ri != Images.end(); ++ri)
    overlapImages.push_back(*ri);

  //sort(overlapImages.begin(), overlapImages.end(), DecreasingArea);
  overlapImages.sort(DecreasingArea);
  for (ReducedImageIterator ri = overlapImages.begin(); ri != overlapImages.end(); )
    {      
      ReducedImage *current = *ri;      
      // at least one arcmin^2 in common with the first one, 
      // should be changed to a common intersection      
      const double overlap = aReference.OverlapArcmin2(*current);
      if (overlap < MinOverlap || !current->IsGoodImage()) 
	{
	  cout <<  " " << current->Name() << " overlap with " << aReference.Name() 
	       << " = " << overlap << " arcmin2" << endl;
	  // this assumes overlapImages is shouldDelete=false;
	  ri = overlapImages.erase(ri);
	}
      else ++ri;
    }
  cout << " Kept " << overlapImages.size() << " images " << endl;
}


Frame UnionFrameWCS(const ReducedImageList& ImList, const ReducedImage* Reference) {

  const ReducedImage* ref = Reference;
  if (!ref) ref = ImList.front();

  GtransfoRef refRaDecToPixels = ref->RaDecToPixels();
  Frame frameRef = ref->UsablePart();

  // initialize
  Frame unionFrame = frameRef;

  for (ReducedImageCIterator it = ImList.begin(); it != ImList.end(); ++it) {
    const ReducedImage *im = *it;
    if (*im == *Reference) continue;
    GtransfoRef imPixelsToRaDec = im->PixelsToRaDec();
    GtransfoRef imToref = GtransfoCompose(refRaDecToPixels, imPixelsToRaDec);
    unionFrame += ApplyTransfo(im->UsablePart(), *imToref, LargeFrame);
  }

  // make an integer frame to avoid resampling1
  unionFrame.xMin = floor(unionFrame.xMin);
  unionFrame.yMin = floor(unionFrame.yMin);
  unionFrame.xMax = ceil(unionFrame.xMax);
  unionFrame.yMax = ceil(unionFrame.yMax);
  cout << " UnionFrame: final frame is " << unionFrame << endl;

  return unionFrame;
}

Frame UnionFrame(const ReducedImageList& ImList, const ReducedImage* Reference) {

  const ReducedImage* ref = Reference;
  if (!ref) ref = ImList.front();

  BaseStarList brefList; LoadForMatch(*ref, brefList);
  Frame unionFrame(ref->UsablePart());

  for (ReducedImageCIterator it = ImList.begin(); it != ImList.end(); ++it) {
    const ReducedImage *im = *it;
    if (*im == *Reference) continue;
    BaseStarList bimList; LoadForMatch(*im, bimList);
    GtransfoRef imToref = FindTransfo(bimList, brefList, *im, *ref);
    Frame frameInRef = ApplyTransfo(im->UsablePart(), *imToref, LargeFrame);
    unionFrame += frameInRef;
  }

  // make an integer frame to avoid resampling1
  unionFrame.xMin = floor(unionFrame.xMin);
  unionFrame.yMin = floor(unionFrame.yMin);
  unionFrame.xMax = ceil(unionFrame.xMax);
  unionFrame.yMax = ceil(unionFrame.yMax);
  cout << " UnionFrame: final frame is " << unionFrame << endl;

  return unionFrame;
}


void MakeUnionRef(const ReducedImageList& ToAlign, const ReducedImage& Reference, const string& unionName, const bool UseWCS) {
  Frame unionFrame;
  if (UseWCS) 
    unionFrame = UnionFrameWCS(ToAlign, &Reference);
  else 
    unionFrame = UnionFrame(ToAlign, &Reference);

  double dx = unionFrame.xMin;
  double dy = unionFrame.yMin;
  GtransfoRef transfoFromRef = new GtransfoLinShift( dx,  dy);
  GtransfoRef transfoToRef   = new GtransfoLinShift(-dx, -dy);
  ImageGtransfo imTransfo(transfoFromRef, transfoToRef, 
			  ApplyTransfo(unionFrame, *transfoToRef, LargeFrame),
			  Reference.Name());  
  TransformedImage transformedRef(unionName, Reference, &imTransfo);
  transformedRef.Execute(ToTransform(Reference));
  FitsHeader headRef(transformedRef.FitsName(), RW);
  headRef.AddOrModKey("PKAUNION", true);
}


string ImageResample(const ReducedImage& Im, const ReducedImage& Ref, const GtransfoRef ImToRef, const GtransfoRef RefToIm) {
  string resampledName = TransformedName(Im.Name(), Ref.Name());
  { 
    ReducedImageRef resampledIm = ReducedImageNew(resampledName);
    if (resampledIm && resampledIm->IsValid() && resampledIm->Execute(ToTransform(Im))) {
      cout << " ImageResample: " << resampledName << " is done\n";
      return resampledName;
    }
  }

  ImageGtransfo *imTransfo;
  if (ImToRef && RefToIm) {
    imTransfo = new ImageGtransfo(RefToIm,
				  ImToRef,
				  Ref.UsablePart(),
				  Ref.Name());
  } else {
    imTransfo = new ImageGtransfo(Ref, Im);
  }
  TransformedImage imResampled(resampledName, Im, imTransfo);
  if (!imResampled.Execute(ToTransform(Im)))
    throw PolokaException(" Failed to produce " + resampledName);
  delete imTransfo;
  return resampledName;
}

string ImageResample(const ReducedImage& Im, const GtransfoRef ImToResampled, const GtransfoRef ResampledToIm) {
  string resampledName = "T_" + Im.Name();
  {
    ReducedImageRef resampledIm = ReducedImageNew(resampledName);
    if (resampledIm && resampledIm->IsValid() && resampledIm->Execute(ToTransform(Im))) {
      cout << " ImageResample: " << resampledName << " is done\n";
      return resampledName;
    }
  }

  Frame origFrame = Im.UsablePart();
  Frame projectedFrame = ApplyTransfo(origFrame,
				      ImToResampled->LinearApproximation(origFrame.Center()));
  // add 5 pixels margin around
  projectedFrame.CutMargin(-5);
  ImageGtransfo *imTransfo = new ImageGtransfo(ResampledToIm,
					       ImToResampled,
					       projectedFrame,
					       Im.Name());
  TransformedImage imResampled(resampledName, Im, imTransfo);
  if (!imResampled.Execute(ToTransform(Im)))
    throw PolokaException(" Failed to produce " + resampledName);
  delete imTransfo;
  return resampledName;
}

string GtransfoName(const DbImage& Ref, const DbImage& Src)
{
  return Src.Dir() + "/gtransfo_to_" + Ref.Name() + ".dat";
}

string ShiftedName(const string &ToShift, const string &Ref)
{
  return "S_" + Ref + "_" + ToShift;
}

static double sign(const double& x) {
  return x>0.? 1. : -1.;
}


string ImageIntegerShift(const ReducedImage& Im, const ReducedImage& Ref, const GtransfoRef ImToRef) {
  string shiftedName = ShiftedName(Im.Name(), Ref.Name());
  {
    ReducedImageRef shiftedIm = ReducedImageNew(shiftedName);
    if (shiftedIm && shiftedIm->IsValid() && shiftedIm->Execute(ToTransform(Im))) {
      cout << " ImageIntegerShift: " << shiftedName << " already produced\n";
      return shiftedName;
    }
  }
  GtransfoRef transfo = ImToRef;
  if (!transfo) transfo = FindTransfo(Im, Ref);
  GtransfoLin approx = transfo->LinearApproximation(Ref.UsablePart().Center());
  GtransfoLin* shift = new GtransfoLin(ceil(approx.Coeff(0,0,0)),
				       ceil(approx.Coeff(0,0,1)),
				       sign(approx.Coeff(1,0,0)),
				       0, 0,
				       sign(approx.Coeff(0,1,1)));
  GtransfoRef refToIm = new GtransfoLin(shift->invert());
  cout << " Transfo from " << Im.Name() << " to " << Ref.Name() << endl;
  cout << *shift;
  GtransfoRef imToRef = shift->Clone();
  delete shift;
  
  ImageGtransfoRef imShift = new ImageGtransfo(refToIm,
					       imToRef,
					       Ref.PhysicalSize(),
					       Ref.Name());
  TransformedImage imShifted(shiftedName, Im, imShift);
  if (!imShifted.Execute(ToTransform(Im)))
    throw PolokaException(" Failed to produce " + shiftedName);

  return shiftedName;
}
