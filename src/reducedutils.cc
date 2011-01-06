#include "reducedutils.h"
#include "starmatch.h"
#include "listmatch.h"
#include "vutils.h" // for DArrayMedian

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
			double &error, const Gtransfo* transfo)
{
  const double distmax(1.);
  BaseStarList *lcur = (BaseStarList *) &CurList;
  BaseStarList *lref = (BaseStarList *) &RefList;
  GtransfoIdentity ident;
  if (!transfo) transfo = &ident;
  StarMatchList *matchlist = ListMatchCollect(*lcur, *lref, transfo, distmax);
  if (matchlist->empty()) 
    {
      cerr << " QuickPhotomRatio() : Error : matchlist empty, returning ratio=1" << endl;
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
  R = PairListMedianPhotomRatio(L);
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
      bool outliers = false;
      for (int k=0; k < count; ++k) if (chi2Vals[k] > chi2Cut) {outliers = true;break;}
      if (!outliers && oldCount == count && chi2Old - chi2 < 1e-3) break;
      chi2Old = chi2;
      oldCount = count;
      if (outliers) niter ++;
      }   while (niter < 10);
  delete [] chi2Vals;
  
  if ( isnan(R))
    {      
      cerr  << "WARNING SlowPhotomRatio : nan value " << endl ;
      return false ;
    }
  cout << " chi2 photom ratio, niter " << chi2/(L.size()-1) << ' ' << niter << endl;
  return (niter<10);
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

  sort(overlapImages.begin(), overlapImages.end(), DecreasingArea);
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


