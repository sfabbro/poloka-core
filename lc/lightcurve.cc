#include <algorithm>  // min_element, copy, for_each
#include <functional> // bind2nd
#include <iterator>   // ostream_iterator
#include <iomanip>    // setw, fixed ...

#include <reducedimage.h>

#include "lightcurve.h"
#include "lcio.h"

// instantiate
template class Fiducial<PhotStar>;

void LightCurve::push_back(const ReducedImage* Rim, const PhotStar *Star)
{
  // look if image is already there
  LightCurve::iterator it = find_if(begin(), end(), bind2nd(mem_fun(&Fiducial<PhotStar>::HasImage), Rim));

  if (it != end())
    {
      cerr << " LightCurve::push_back() : Warning: image "
	   << Rim->Name() << " is already in the LightCurve " << endl;
      return;
    }

  Fiducial<PhotStar> *fidStar = new Fiducial<PhotStar>(Star, Rim);

  push_back(fidStar);
}

void LightCurve::push_back(const ReducedImage* Rim)
{
  // look if image is already there
  LightCurve::iterator it = find_if(begin(), end(), bind2nd(mem_fun(&Fiducial<PhotStar>::HasImage), Rim));

  if (it != end())
    {
      cerr << " LightCurve::push_back() : Warning: image "
	   << Rim->Name() << " is already in the LightCurve " << endl;
      return;
    }

  Fiducial<PhotStar> *fidStar = new Fiducial<PhotStar>(Rim);

  push_back(fidStar);
}

void LightCurve::write_short(ostream& Stream) const
{
  ios::fmtflags oldflags = Stream.flags();

  if (front()->Image()) Stream << "# jd \n";
  Stream << "# flux \n"  
         << "# eflux \n";
  if (front()->Image()) Stream << "# image \n";

  Stream << setiosflags(ios::fixed);
  for (LightCurve::const_iterator it = begin(); it != end(); ++it)
    {      
      const Fiducial<PhotStar> *fs = *it;
      if (fs->Image()) Stream << setw(14) << setprecision(2) << fs->Image()->JulianDate();
      else Stream << 99999.99;
      Stream << setw(15) << setprecision(3) << fs->flux
	     << setw(15) << setprecision(3) << sqrt(fs->varflux);
      if (fs->Image()) Stream << "  " << fs->Image()->Name();
      else Stream << " none ";
      Stream << endl;
    }

  Stream.flags(oldflags);
}


ostream& operator << (ostream& Stream, const CountedRef<Fiducial<PhotStar> > &Star)
{
  Stream << *Star;
  return Stream;
}

ostream& operator << (ostream& Stream, const LightCurve& Lc)
{
  copy(Lc.begin(), Lc.end(), ostream_iterator<CountedRef<Fiducial<PhotStar> > >(Stream, "\n"));
  return Stream;
}



LightCurveList::LightCurveList(istream& LcFileStream) 
{

  lc_read(LcFileStream, Objects, Images);
  
  RefImage = Objects.front()->Image();

  // fill up the light curve
  for (RefStarCIterator it = Objects.begin(); it != Objects.end(); ++it)
    {
      LightCurve lc(*it);
      // foreach object, link the list of images
      for (ReducedImageCIterator im=Images.begin(); im != Images.end(); ++im) {
	PhotStar *fidPhot = new PhotStar(BaseStar((*it)->x, (*it)->y, 0.));
	lc.push_back(*im, fidPhot); // add one image and one PhotStar
      }
      push_back(lc);
    }  
  
  

  cout << " LightCurveList::LightCurveList() : " << size() 
       << " objects " << Images.size() << " images \n";
}

ostream& operator << (ostream& Stream, const LightCurveList& Fiducials)
{
  copy(Fiducials.begin(), Fiducials.end(), ostream_iterator<LightCurve>(Stream));
  return Stream;
}
