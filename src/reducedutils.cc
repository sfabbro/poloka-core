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

static double sqr(const double &x) {return x*x;}

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
	  const double weight = 1./(sqr(ratio*(cur->EFlux()/cur->flux + ref->EFlux()/ref->flux)));
	  sumwt += weight;
	  sum += weight * ratio;
	}
    }
  error = 1/sqrt(sumwt);
  return sum / sumwt;
}



double MedianPhotomRatio(const StarMatchList *matchlist)
{
  if (matchlist->empty()) 
    {
      cerr << " MedianPhotomRatio() : Error: matchlist empty, returning ratio=-99" << endl;
      return -99.;
    }


  double *ratios = new double[matchlist->size()];
  int count = 0;
  for (StarMatchCIterator i = matchlist->begin(); i != matchlist->end(); ++i)
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
  double pr= MedianPhotomRatio(matchlist);
  delete matchlist;
  return pr;
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


