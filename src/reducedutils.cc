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

static bool SM_DecFluxMax(const StarMatch &S1, const StarMatch &S2)
{
  const SEStar *p1 = (const SEStar *) (const BaseStar *) S1.s1;
  const SEStar *p2 = (const SEStar *) (const BaseStar *) S2.s1;
  return  DecFluxMax(p1,p2);
}

/************************ Photom ratio estimators *******************/

/*  notes:

Measuring a flux ratio is not that straightforward. The "obvious"
least-squares approach is for example facing serious problems.
Imagine you want to fit y=ax on a sample of x's and y's, i.e. estimate the
ratio a. If you minimize sum_i(y_i-a * x_i)^2/sy_i^2, the estimate
of a will be unbiased only if there are no uncertainties on x. In the opposite
case, there is a bias : <a_min> = <a_true> * sum(x_i^2)/ sum(x_i^2+sx_i^2).
There is no single way to get around this problem. Getting back to photometric
ratios, there is usually no problem about optimality, because the 
statistical accuracy is usually far too good, compared to systematics
(e.g. non-uniformity of response due to PSF and/or flatfielding).
  QuickPhotomRatio computes a weighted average of the ratios. This
is less biased than straight LS, but still biased because <1/x> != 1/<x>
for symetric pdf's. One improvement is to compute both ratios (r1 and r2);
if the estimator is unbaised, then <r1*r2> = 1; one way to compensate
for a bias of r1 and r2 is to return the geometric mean sqrt(r1/r2);
This is not coded (yet).

  One estimator with a very small bias, and almost optimal is
to use the functional relationship framework. Namely, you use least-squares:

chi2 = sum (y_i-a xp_i)^2/sy_i^2 + (x_i-xp_i)^2/sx_i^2

where a and xp_i are parameters. When you solve the normal equations,
you find that it is equivalent to minimize:

chi2 = sum (y_i - a x_i)^2/(sy_i^2+(a x_i)^2)

The minimization is tedious (non linearities are important), but
doable. This is done in SlowPhotomRatio()


compute

*/

 

static double sq(const double &x) {return x*x;}

double QuickPhotomRatio(const ReducedImage &CurImage, const ReducedImage &RefImage, 
			double &error, const Gtransfo* cur2ref)
{
  if (CurImage == RefImage) return 1.;

  SEStarList curList(CurImage.CatalogName());
  SEStarList refList(RefImage.CatalogName());

  curList.FluxSort();
  curList.CutTail(200);
  refList.FluxSort();
  refList.CutTail(200);

  return QuickPhotomRatio(curList, refList, error, cur2ref);
}

double QuickPhotomRatio(const SEStarList &CurList, const SEStarList &RefList, 
			double &error, const GtransfoRef transfo)
{
  const double distmax(1.);
  BaseStarList *lcur = (BaseStarList *) &CurList;
  BaseStarList *lref = (BaseStarList *) &RefList;
  GtransfoIdentity ident;
  const Gtransfo* mytransfo = transfo;
  if (!transfo) mytransfo = &ident;
  StarMatchList *matchlist = ListMatchCollect(*lcur, *lref, mytransfo, distmax);
  if (matchlist->empty()) 
    {
      cerr << " QuickPhotomRatio() : Error : matchlist empty, returning ratio=1" << endl;
      error = 0;
      return 1;
    }
  matchlist->sort(&SM_DecFluxMax);
  int count = 0 ; 
  StarMatchIterator it = matchlist->begin();
  int nmax = 200;
  while ((it != matchlist->end()) && ( count < nmax)) {++count; ++it;}
  matchlist->erase(it, matchlist->end()) ;
  double sumwt=0,sum=0;
  for (it = matchlist->begin(); it != matchlist->end(); ++it) 
    {      
      const SEStar *cur = it->s1;
      const SEStar *ref = it->s2;
      if (ref->flux>0.0 && cur->flux>0.0 && cur->EFlux()>0.0 && ref->EFlux()>0.0)
	{
	  const double ratio = cur->flux / ref->flux;
	  const double weight = 1./(sq(ratio*(cur->EFlux()/cur->flux + ref->EFlux()/ref->flux)));
	  sumwt += weight;
	  sum += weight * ratio;
	}
    }
  error = 1/sqrt(sumwt);
  return sum / sumwt;
}



double MedianPhotomRatio(const StarMatchList &MatchList)
{
  if (MatchList.empty()) 
    {
      cerr << " MedianPhotomRatio() : Error: matchlist empty, returning ratio=-99" << endl;
      return -99.;
    }
  /* This is unweighted. If we weight it, we might consider 
     averaging between the direct and reversed photom ratio 
  */
  double *ratios = new double[MatchList.size()];
  int count = 0;
  for (StarMatchCIterator i = MatchList.begin(); i != MatchList.end(); ++i)
    if(i->s2->flux>0.)
      ratios[count++] = i->s1->flux/i->s2->flux;
  double med = DArrayMedian(ratios,count);
  delete [] ratios;
  return med;
}

double MedianPhotomRatio(const BaseStarList &CurList, const BaseStarList &RefList, const Gtransfo *Transfo)
{
  StarMatchList *matchlist;
  if (!Transfo) matchlist = ListMatchCollect(CurList, RefList, 1.);
  else matchlist = ListMatchCollect(CurList, RefList, Transfo, 1.);
  /* it is useless to compute both since one is exactly the 
     inverse of the other : */
#ifdef STORAGE
  double pr1 = MedianPhotomRatio(*matchlist);
  matchlist->Swap();
  double pr2 = MedianPhotomRatio(*matchlist);
  delete matchlist;
  return sqrt(pr1/pr2);
#endif
  double pr = MedianPhotomRatio(*matchlist);;
  delete matchlist;
  return pr;
}


static double PairListMedianPhotomRatio(const FluxPairList &L)
{
  double *ratio = new double[L.size()];
  int count = 0;
  for (FluxPairList::const_iterator i = L.begin(); i != L.end(); ++i)
    {
      ratio[count++] = i->f1/i->f2;
    }
  double pr = DArrayMedian(ratio,count);
  delete [] ratio;
  return pr;
}

//! see notes above if you wonder what is done here.
bool SlowPhotomRatio(const FluxPairList &L, const double NSigChi2Cut, double &R, double &Var)
{
  // compute a such that f1=a*f2 on average.
  int niter = 0;
  double chi2;
  double chi2Old = 1e30;
  int oldCount = L.size();
  double *chi2Vals = new double[oldCount];
  double chi2Cut = 1e30;
  const int nitermax = 50;
  R = PairListMedianPhotomRatio(L);
  // cout << " initial ratio = " << R << " with " << L.size() << " stars" << endl;
  do
    {
      double num = 0;
      double deno = 0;
      chi2 = 0;
      int count = 0;
      for (FluxPairList::const_iterator i = L.begin(); i != L.end(); ++i)
	{
	  double y = i->f1;
	  double x = i->f2;
	  double sx2 = sq(i->sig2);
	  double sy2 = sq(i->sig1);
	  double d = (sy2+R*R*sx2);
	  double xp = (R*y*sx2+x*sy2)/d;
	  double thisChi2 = sq(y-R*x)/d;
	  if (thisChi2 > chi2Cut) continue;
	  chi2Vals[count++] = thisChi2;
	  chi2 += thisChi2;
	  num += xp*(y-R*xp)/sy2;
	  deno += (sq(xp)-sx2*sq(y-2.*R*xp)/d)/sy2;

	}
      // no protection for dividing by zero
      if ( fabs(deno) < 1.e-20 )
	{
	  cerr << "WARNING SlowPhotomRatio : dividing by zero" << endl ; 
	  return false ;
	}
      R += num/deno;
      Var = 1/deno;
      double chi2Mean, chi2Med, chi2Sig;
      Dmean_median_sigma(chi2Vals, count, chi2Mean, chi2Med, chi2Sig);
      chi2Cut = chi2Med+NSigChi2Cut*chi2Sig;
      if (chi2Cut == 0)
	{
	  cout << " INFO : SlowPhotomRatio finds a null chi2 cut !! seems we got twice the same list " << endl;
	  chi2Cut = 1;
	}
      bool outliers = false;
      for (int k=0; k < count; ++k) if (chi2Vals[k] > chi2Cut) {outliers = true;break;}
      if (!outliers && oldCount == count && chi2Old - chi2 < 1e-3) break;
      chi2Old = chi2;
      oldCount = count;
      if (outliers) niter++;
      }   while (niter < nitermax);
  delete [] chi2Vals;

  if ( isnan(R))
    {      
      cerr  << " WARNING SlowPhotomRatio : nan value " << endl ;
      return false ;
    }

  // cout << " chi2 photom ratio, niter " << chi2/(L.size()-1) << ' ' << niter << endl;
  return (niter<nitermax);
}

void LoadForMatch(const ReducedImage& Im, BaseStarList& BList, const double& MinSN) {
  
  AperSEStarList sList;
  if (Im.HasAperCatalog())
    sList.read(Im.AperCatalogName());
  else if (Im.HasCatalog())
    sList.read(Im.CatalogName());
  
  for (AperSEStarIterator it=sList.begin(); it != sList.end(); ) {
    AperSEStar *pstar = *it ;
    // not in bad zone, keep satur, only with better 0.3 pixels accuracy
    if (pstar->FlagBad() == 0 && pstar->Flag() < 8 && 
	pstar->flux > MinSN * pstar->EFlux())
      //pstar->vx + pstar->vy < 0.1)
      ++it;
    else 
      it = sList.erase(it);
  }  
  BList = *AperSE2Base(&sList);
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
  size_t nmin = size_t(min(Src.size(), Dest.size()) * MinMatch);
  cout << " ListMatchCheck: min stars to match = " << nmin << endl;
  size_t ngood = StarMatchCheck(match);
  delete match;
  return (ngood > nmin);
}

GtransfoRef FindWCSTransfo(const ReducedImage& Src, const ReducedImage& Dest) {

  GtransfoRef src2dest;
  FitsHeader hsrc(Src.FitsName());
  if (!HasLinWCS(hsrc)) return src2dest;
  FitsHeader hdest(Dest.FitsName());
  if (!HasLinWCS(hdest)) return src2dest;
  GtransfoRef wsrc = WCSFromHeader(hsrc);
  GtransfoRef wdest = WCSFromHeader(hdest);
  
  if (wsrc && wdest) {
    Frame fsrc = Src.UsablePart();
    Frame fdest = Dest.UsablePart();
    src2dest = GtransfoCompose(wdest->InverseTransfo(0.01, fdest), wsrc);
  }
  return src2dest;
}

GtransfoRef FindTransfo(const ReducedImage& Src, const ReducedImage& Dest, int MaxOrder) {

  if (Src == Dest) {
    cout << " FindTransfo: same image: " << Src.Name() << " returning identity\n";
    return new GtransfoIdentity();
  }
    
  BaseStarList bsrcList;  LoadForMatch(Src, bsrcList);
  BaseStarList bdestList; LoadForMatch(Dest, bdestList);

  return FindTransfo(bsrcList, bdestList, Src, Dest, MaxOrder);  
}

GtransfoRef FindTransfo(const BaseStarList& SrcList, const BaseStarList& DestList,
			const ReducedImage& Src, const ReducedImage& Dest,
			int MaxOrder) {

  cout << " FindTransfo: " 
       << Src.Name() << " (" << SrcList.size() << " stars) to " 
       << Dest.Name() << " (" << DestList.size() << " stars)\n";
  cout << " FindTransfo: trying with WCS's...\n";
  GtransfoRef transfo = FindWCSTransfo(Src, Dest);
  MatchConditions cond(DefaultDatacards());
  cond.SizeRatio = Src.PixelSize() / Dest.PixelSize();
  cond.DeltaSizeRatio = 0.1 * cond.SizeRatio;

  if (transfo) {
    if (ListMatchCheck(SrcList, DestList, transfo, 2, cond.MinMatchRatio)) {
      cout << " FindTransfo: refining WCS transfo\n";
      transfo = ListMatchRefine(SrcList, DestList, transfo, MaxOrder);
    } else {
      cout << " FindTransfo: WCS not good enough\n";
      transfo = GtransfoRef();
    }
  }

  if (!transfo) {
    cout << " FindTransfo: trying with combinatorial match\n";
    transfo = ListMatchCombinatorial(SrcList, DestList, cond);
    if (transfo) {
      cout << " FindTransfo: refining combinatorial match\n";
      transfo = ListMatchRefine(SrcList, DestList, transfo);
    } else {
      cout << " FindTransfo: no match found\n";
      transfo = GtransfoRef();
    }
  }
  return transfo;
}

bool PhotomRatio(const DbImage &Im, const DbImage &Ref, double& Ratio, double &Error, GtransfoRef &Im2Ref)
{
  const size_t nmax = 500;
  if (Im == Ref) { 
    Error = 0.;
    Ratio = 1.;
    return true;
  }
  string imListName = Im.ImageCatalogName();
  string refListName = Ref.ImageCatalogName();

  if (!FileExists(imListName) || !FileExists(refListName)) {
   imListName = Im.ImageCatalogName(SExtractor);
   refListName = Ref.ImageCatalogName(SExtractor);
   if (!FileExists(imListName) || !FileExists(refListName)) {
     cerr << " Error: no similar catalogues available for "
	  << Im.Name() << " or " << Ref.Name() << endl;
     return false;
   }
  }

  SEStarList imList(imListName);
  SEStarList refList(refListName);
  BaseStarList *bimList  = SE2Base(&imList);
  BaseStarList *brefList = SE2Base(&refList); 
  // hack to detect a transformed image
  if (!Im2Ref)
    Im2Ref = ListMatch(*bimList, *brefList);
  if (!Im2Ref) {
     cerr << " Error: bad match for PhotomRatio between "
	  << Im.Name() << " and " << Ref.Name() << endl;
     return false;
  }

  StarMatchList *matchList = ListMatchCollect(*bimList, *brefList, Im2Ref, 1.);
  matchList->sort(&SM_DecFluxMax);
  FluxPairList fpl;
  for (StarMatchCIterator it = matchList->begin(); it != matchList->end() && fpl.size() < nmax; ++it)  {
    const SEStar* s = (const SEStar *)(it->s1);
    const SEStar* r = (const SEStar *)(it->s2);
    if (!r->IsBad() && !r->IsSaturated() && r->flux > 0 && r->EFlux() > 0 &&
	!s->IsBad() && !s->IsSaturated() && s->flux > 0 && s->EFlux() > 0)
      fpl.push_back(FluxPair(s->flux, s->EFlux(),r->flux, r->EFlux()));
  }
  delete matchList;
  if (fpl.empty()) {
    cerr << " Error: empty list in PhotomRatio between "
	 << Im.Name() << " and " << Ref.Name() << endl;
    return false;
  }
  double var;
  bool ok = SlowPhotomRatio(fpl, 5, Ratio, var);
  Error = sqrt(var);
  if (!ok) {
    cerr << " Error: could not find photometric ratio for " 
	 << Im.Name() << " vs. " << Ref.Name() << endl;
    ok = false;
  }
  return ok;
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

void MakeUnionRef(const ReducedImageList& ToAlign, const ReducedImage& Reference, const string& unionName) {
  Frame unionFrame = UnionFrame(ToAlign, &Reference);
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


string ImageResample(const ReducedImage& Im, const ReducedImage& Ref) {
  string resampledName = TransformedName(Im.Name(), Ref.Name());
  { 
    ReducedImageRef resampledIm = ReducedImageNew(resampledName);
    if (resampledIm && resampledIm->IsValid() && resampledIm->Execute(ToTransform(Im))) {
      cout << " ImageResample: " << resampledName << " already produced\n";
      return resampledName;
    }
  }
  GtransfoRef imToRef = FindTransfo(Im, Ref);
  GtransfoRef refToIm = imToRef->InverseTransfo(0.001, Ref.UsablePart());
  cout << " Transfo from " << Im.Name() << " to " << Ref.Name() << endl;
  cout << *imToRef;
  ImageGtransfo *imTransfo = new ImageGtransfo(refToIm,
					       imToRef,
					       Ref.PhysicalSize(),
					       Ref.Name());
  TransformedImage imResampled(resampledName, Im, imTransfo);
  if (!imResampled.Execute(ToTransform(Im)))
    throw PolokaException(" Failed to produce " + resampledName);
  delete imTransfo;
  return resampledName;
}

static string ShiftedName(const string &ToShift, const string &Ref)
{
  return "S_" + Ref + ToShift;
}

static double sign(const double& x) {
  return x>0.? 1. : -1.;
}


string ImageIntegerShift(const ReducedImage& Im, const ReducedImage& Ref) {
  string shiftedName = ShiftedName(Im.Name(), Ref.Name());
  {
    ReducedImageRef shiftedIm = ReducedImageNew(shiftedName);
    if (shiftedIm && shiftedIm->IsValid() && shiftedIm->Execute(ToTransform(Im))) {
      cout << " ImageIntegerShift: " << shiftedName << " already produced\n";
      return shiftedName;
    }
  }
  GtransfoRef transfo = FindTransfo(Im, Ref);
  GtransfoLin lintransfo = transfo->LinearApproximation(Ref.UsablePart().Center());
  GtransfoLin *imToRef = new GtransfoLin(ceil(lintransfo.Coeff(0,0,0)),
					 ceil(lintransfo.Coeff(0,0,1)),
					 sign(lintransfo.Coeff(1,0,0)),
					 0, 0,
					 sign(lintransfo.Coeff(0,1,1)));

  GtransfoRef refToIm = new GtransfoLin(imToRef->invert());
  cout << " Transfo from " << Im.Name() << " to " << Ref.Name() << endl;
  cout << *imToRef;

  ImageGtransfo *imShift = new ImageGtransfo(refToIm,
					     imToRef,
					     Ref.PhysicalSize(),
					     Ref.Name());
  TransformedImage imShifted(shiftedName, Im, imShift);
  if (!imShifted.Execute(ToTransform(Im)))
    throw PolokaException(" Failed to produce " + shiftedName);
  delete imShift;
  delete imToRef;
  return shiftedName;
}
