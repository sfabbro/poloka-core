#include <iomanip>
#include "reducedutils.h"
#include "fileutils.h"
#include "imagesum.h"
#include "starmatch.h"
#include "imagematch.h"
#include "listmatch.h"

#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif

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

  double ratio =  QuickPhotomRatio(curList, refList, error, cur2ref);
  cout << " ratio  "<< CurImage.Name() << " / " << RefImage.Name() 
       << " = " << ratio << " +/- " << error << endl; 

  return ratio;
}

double QuickPhotomRatio(const SEStarList &CurList, const SEStarList &RefList, 
			double &error, const Gtransfo* transfo)
{
  double distmax = 1.;
  BaseStarList *lcur = (BaseStarList *) &CurList;
  BaseStarList *lref = (BaseStarList *) &RefList;
  GtransfoIdentity ident;
  if (!transfo) transfo = &ident;
  StarMatchList *matchlist = ListMatchCollect(*lcur, *lref, transfo, distmax);
  if (matchlist->size() == 0 ) 
    {
      cerr << " FAILURE QuickPhotomRatio : matchlist empty" << endl; 
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
	  double ratio = cur->flux / ref->flux;
	  double weight = 1./(sqr(ratio*(cur->EFlux()/cur->flux + ref->EFlux()/ref->flux)));
	  sumwt += weight;
	  sum += weight * ratio;
	}
    }
  error = 1/sqrt(sumwt);
  return sum / sumwt;
}

NightSet::NightSet(const ReducedImage &OneImage)
  :ReducedImageList(true)
{
  push_back(OneImage.Clone());
}

NightSet::NightSet(const NightSet &Other) : ReducedImageList(Other) {};

bool NightSet::IsFromSameSet(const ReducedImage &Another) const
{  
  if (find(begin(), end(), &Another) != end()) return false;
  string curdate = front()->Date();
  string othdate = Another.Date();
  return ((front()->SameChipFilterInst(Another,false)) && (othdate==curdate));
}

string NightSet::GenericName() const
{
  //re-string the date
  string datestr = front()->Date();
  int dd,mm,yy;
  sscanf(datestr.c_str(),"%d/%d/%d",&dd,&mm,&yy);
  char datechr[64];
  sprintf(datechr,"%d%d%d",yy,mm,dd);

#ifdef ALOTOFIMAGES
  //remove slashes, commas and blanks from Target field
  string objestr = front()->Target();
  string badchars = "/ ,;-";
  for (unsigned i=0; i< badchars.length(); ++i) 
    {
      int s = objestr.find(badchars[i]);
      if (s>0) objestr.erase(s,objestr.length());
    }
  string genname = string(datechr)+"_"+objestr+"_"+front()->Telescope()+"_"+front()->Band();
#endif

  string genname = front()->Telescope() + string(datechr)+ "_" + front()->Band();
  return genname;
}

double NightSet::Seeing() const
{
  unsigned n=size();
  double seeing = 0;
  for (unsigned i=0;i<n;++i) seeing += (*this)[i]->Seeing();
  return seeing/n;
}

double NightSet::Exposure() const
{
  double expo = 0;
  for (unsigned i=0;i<size();++i) expo += (*this)[i]->Exposure();
  return expo;
}

double NightSet::SignalToNoise23() const
{
  double signal = 0;
  double noise2 = 0;
  for (unsigned i=0;i<size();++i) 
    {
      signal += pow(10., -(23.- (*this)[i]->ZeroPoint())/2.5);
      noise2 += M_PI*sqr(2.3548 * (*this)[i]->Seeing() * (*this)[i]->SigmaBack());
      
    }
  return signal/sqrt(noise2);
}

void NightSet::dump(ostream& Stream) const
{
  Stream << endl;
  Stream << "**************NightSet composition*************************" << endl ;
  Stream << "  Date: " << front()->Date() << "  Target: " << front()->Target() << endl
	 << "  Telescope: " << front()->Telescope()  << "  Instrument: " << front()->Instrument()
	 << "  Filter: " << front()->Band() << endl;
  Stream << "-----------------------------------------------------------" << endl;
  Stream << setw(17) << setiosflags(ios::left)<< "Image" << resetiosflags(ios::left);
  Stream << setw(10) << "Exposure" << setw(10) 
  	 << "Seeing" << setw(8) << "Sky" << setw(12) << "S/N(23)" << endl;
  Stream << "-----------------------------------------------------------" << endl;
  Stream << setiosflags(ios::fixed);

  for (ReducedImageCIterator ri = begin() ; ri!=end(); ++ri)
    {
      const ReducedImage *current = *ri;
      Stream << setw(17) << setiosflags(ios::left)<< current->Name() 
	     << resetiosflags(ios::left);
      Stream << setw(10) << setprecision(0) << current->Exposure()
	     << setw(10) << setprecision(3) << current->Seeing()*2.3548*current->PixelSize()
	     << setw(10) << setprecision(0) << current->BackLevel()
	     << setw(10) << setprecision(2) << current->SignalToNoise23() << endl;
    }
  Stream << resetiosflags(ios::fixed);
}


// to arrange a list night by night
void ArrangeByNight(const ReducedImageList &Images, vector<NightSet> &AllNights)
{
  cout << " Arranging images by night " << endl;
  ReducedImageList imlist(false);
  for (ReducedImageCIterator ri = Images.begin() ; ri!=Images.end(); ++ri) 
    imlist.push_back(*ri);

  // image by image empty the imageslist
  while (!imlist.empty())
    {
      // get the first one to make a new NightSet
      ReducedImage *first = imlist.front();
      NightSet thisnight(*first);      
      imlist.erase(imlist.begin());

      // check all others, copy and add entry if from same set delete otherwise
      for (ReducedImageIterator ri = imlist.begin() ; ri!=imlist.end(); )
	{
	  ReducedImage *current = *ri;
	  if (thisnight.IsFromSameSet(*current)) 
	    {
	      thisnight.push_back(current->Clone()); 
	      ri = imlist.erase(ri);
	    }
	  else ++ri;
	}
      AllNights.push_back(thisnight);
    }
  cout << " We are dealing with " << AllNights.size() << " night sets" << endl;
}

void SelectBestImages(ReducedImageList &Images)
{
  cout << " Select best images " << endl;
  // first select images with smallest pixel size (5% tolerance)
  sort(Images.begin(), Images.end(), IncreasingResolution);
  double maxpixsize = Images.front()->PixelSize()*1.05;
  ReducedImageIterator ri = Images.begin();
  while (ri != Images.end() && (*ri)->PixelSize()< maxpixsize) ++ri;
  ri = Images.erase(ri, Images.end());
									
  // now select images of best S/N and best seeing of a given instrument
  double bestseeing = 10;
  ReducedImageList bestImages = Images;

  while (!bestImages.empty())
    {
      // compute seeing for images of a same telescope
      double seeing = 0;
      ReducedImageList curList(false);
      ReducedImageIterator it = bestImages.begin();
      ReducedImage *first = *it;
      while (it != bestImages.end())
	{
	  ReducedImage *current = *it;
	  if (current->SameChipFilterInst(*first,false)) 
	    {
	      double curseeing = current->Seeing();
	      seeing += curseeing;
	      curList.push_back(current);
	      it = bestImages.erase(it);
	    }
	  else ++it;
	}
      seeing /= curList.size();
      if (seeing < bestseeing)
	{
	  bestseeing = seeing;
	  Images = curList;
	}
      curList.clear();
    }

  double seeingRatio = 1.1;
  if (getenv("SEEINGSPAN")) seeingRatio = atof(getenv("SEEINGSPAN"));
  sort(Images.begin(), Images.end(), IncreasingSeeing);
  double maxseeing = Images.front()->Seeing()*seeingRatio;
  ri = Images.begin();
  while (ri != Images.end() && (*ri)->Seeing()< maxseeing) ++ri;
  ri = Images.erase(ri, Images.end());
  for (ri = Images.begin(); ri != Images.end(); ++ri) 
    cout << " " << (*ri)->Name() << " seeing " << (*ri)->Seeing() << endl;
}

const ReducedImage *BuildReference(const string& RefName, const ReducedImageList &Images, const bool ToAlign)
{  
  ReducedImage red(RefName);
  if (red.ActuallyReduced()) return red.Clone();
  ReducedImageList toStack(false);
  for (ReducedImageCIterator ri = Images.begin() ; ri != Images.end(); ++ri) 
    toStack.push_back(*ri);
  SelectBestImages(toStack);
  ImageSum *ref = new ImageSum(RefName, toStack);
  int toDo = DoFits | DoCatalog;
  if (ToAlign)
    {
      delete ref;
      ref = ImagesAlignAndSum(toStack, *toStack.front(), RefName, toDo);
    }
  else ref->Execute(toDo);
  return ref;
}

// to find a good geometric reference
const ReducedImage* BestResolutionReference(const ReducedImageList &ImList)
{
  cout << " Now looking for the best resolution reference " << endl;  
  const ReducedImage *geomRef = ImList.front();
  //loop over the list to find it
  for (unsigned int i=0; i < ImList.size() ; ++i)
    {      
      const ReducedImage *current = ImList[i];
      double pixratio = current->PixelSize() / geomRef->PixelSize();
      // if pixel size is 10% smaller, keep it. Otherwise, take best seeing
      if (pixratio < 0.9) geomRef = current;
      else if (current->Seeing() <= geomRef->Seeing() && pixratio <1.01) 
	geomRef = current;
    }
  cout << " Found " << geomRef->Name() << ": seeing " << geomRef->Seeing() 
       << "  pixel " << geomRef->PixelSize() << endl;
  
  //do not return a copy of it, will be deleted when ImList is
  return geomRef;
}

// to get rid of images non overlapping in a list of images
void FilterWithOverlap(const ReducedImageList &Images, const ReducedImage &aReference, 
		       ReducedImageList &overlapImages, const double &MinOverlap)
{
  for (ReducedImageCIterator ri = Images.begin(); ri != Images.end(); ++ri)
    overlapImages.push_back(*ri);

  sort(overlapImages.begin(), overlapImages.end(), DecreasingArea);
  for (ReducedImageIterator ri = overlapImages.begin(); ri != overlapImages.end(); )
    {      
      ReducedImage *current = *ri;      
      // at least one arcmin^2 in common with the first one, 
      // should be changed to a common intersection      
      double overlap = aReference.OverlapArcmin2(*current);
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


