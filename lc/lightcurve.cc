#include <algorithm>  // min_element, copy, for_each
#include <functional> // bind2nd
#include <iterator>   // ostream_iterator
#include <iomanip>    // setw, fixed ...

#include <fitsimage.h>
#include <reducedimage.h>

#include "lightcurve.h"
#include "lcio.h"

#define FNAME

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

  double elixir_zp = computeElixirZeroPoint();

  ios::fmtflags oldflags = Stream.flags();

  if (front()->Image()) Stream << "# jd : \n";
  Stream << "# flux : \n"  
         << "# eflux : \n"
	 << "# mag : using elixir zero point = " << elixir_zp << "\n"
	 << "# emag_minus : useful for drawing\n"
    	 << "# emag_plus : useful for drawing\n";
  
  
  if (front()->Image()) Stream << "# image : \n";
  Stream << "# end \n";

  Stream << setiosflags(ios::fixed);
  for (LightCurve::const_iterator it = begin(); it != end(); ++it)
    {      
      const Fiducial<PhotStar> *fs = *it;
      if (fs->Image()) Stream << setw(14) << setprecision(2) << fs->Image()->JulianDate();
      else Stream << 99999.99;
      Stream << setw(15) << setprecision(3) << fs->flux
	     << setw(15) << setprecision(3) << sqrt(fs->varflux);
      if(fs->flux < 1.e-12) {
	Stream << setw(15) << setprecision(3) << 99
	       << setw(15) << setprecision(3) << 0
	       << setw(15) << setprecision(3) << 0;
      }else{
	Stream << setw(15) << setprecision(3) << -2.5*log10(fs->flux)+elixir_zp
	       << setw(15) << setprecision(3) << 2.5*log10(1.-sqrt(fs->varflux)/fs->flux)
	       << setw(15) << setprecision(3) << 2.5*log10(1.+sqrt(fs->varflux)/fs->flux);
      }
      
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
#ifdef FNAME
  cout << " > LightCurveList::LightCurveList(istream& LcFileStream)" << endl;
#endif
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

double LightCurve::computeElixirZeroPoint() const {
  // > COMMENT   Formula for Photometry, based on keywords given in this header:
  // > COMMENT   m = -2.5*log(DN) + 2.5*log(EXPTIME)
  // > COMMENT   M = m + PHOT_C + PHOT_K*(AIRMASS - 1) + PHOT_X*(PHOT_C1 - PHOT_C2)
  
  // we do not take into account color terms !!!
  FitsHeader refhead(Ref->Image()->FitsName());
  double expo =  refhead.KeyVal("TOADEXPO");
  double PHOT_C =  refhead.KeyVal("PHOT_C");
  double PHOT_K =  refhead.KeyVal("PHOT_K");
  double AIRMASS =  refhead.KeyVal("AIRMASS");
  
  return 2.5*log10(expo) + PHOT_C + PHOT_K*(AIRMASS-1.);
}
