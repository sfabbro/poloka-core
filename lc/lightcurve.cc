#include <algorithm>  // min_element, copy, for_each
#include <functional> // bind2nd
#include <iterator>   // ostream_iterator
#include <iomanip>    // setw, fixed ...

#include <fitsimage.h>
#include <reducedimage.h>

#include "lightcurve.h"
#include "lcio.h"

// for io
#include "lightcurvepoint.h"
#include "lightcurvepoint_dict.h"
#include "objio.h"
#include "typemgr.h"


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
	Stream << setw(15) << setprecision(3) << -2.5*log10(fs->flux)+elixir_zp;
	if((1.-sqrt(fs->varflux)/fs->flux)<0)
	  Stream  << setw(15) << setprecision(3) << -99;
	else
	  Stream  << setw(15) << setprecision(3) << 2.5*log10(1.-sqrt(fs->varflux)/fs->flux);
	Stream  << setw(15) << setprecision(3) << 2.5*log10(1.+sqrt(fs->varflux)/fs->flux);
      }
      
      if (fs->Image()) Stream << "  " << fs->Image()->Name();
      else Stream << " none ";
      Stream << endl;
    }

  Stream.flags(oldflags);
}

void LightCurve::write_xml(const string &filename) const
{
#ifdef FNAME
  cout << " > LightCurve::write_xml" << endl;
#endif

  double elixir_zp = computeElixirZeroPoint();
  
  // fill a list of LightCurvePoint
  std::list< CountedRef<LightCurvePoint> > lcpoints;  
  for (LightCurve::const_iterator it = begin(); it != end(); ++it)
    {      
      const Fiducial<PhotStar> *fs = *it;
      CountedRef<LightCurvePoint> lcp = new LightCurvePoint();
      lcp->julianday = fs->Image()->JulianDate();
      lcp->flux = fs->flux;
      lcp->eflux = sqrt(fs->varflux);
      if(fs->flux < 1.e-12) {
	lcp->mag = 99;
	lcp->emag_minus = 0;
	lcp->emag_plus = 0;
      }else{
	lcp->mag = -2.5*log10(fs->flux)+elixir_zp;
	if((1.-sqrt(fs->varflux)/fs->flux)<0)
	  lcp->emag_minus = -99;
	else
	  lcp->emag_minus = 2.5*log10(1.-sqrt(fs->varflux)/fs->flux);
	lcp->emag_plus = 2.5*log10(1.+sqrt(fs->varflux)/fs->flux);
      }
      lcpoints.push_back(lcp);
    }
  // now write this list in a file
  obj_output<xmlstream> oo(filename);
  oo << lcpoints;
  oo.close();
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

#define USE_SKMAGATT

double LightCurve::computeElixirZeroPoint() const {
  // > COMMENT   Formula for Photometry, based on keywords given in this header:
  // > COMMENT   m = -2.5*log(DN) + 2.5*log(EXPTIME)
  // > COMMENT   M = m + PHOT_C + PHOT_K*(AIRMASS - 1) + PHOT_X*(PHOT_C1 - PHOT_C2)
  
  
  string photometric_image_fitsname = Ref->Image()->FitsName(); // default
  double photomratio = 1; // flux(photometric_image)/flux(reference_image=Ref->Image())
  
#ifdef USE_SKMAGATT

  // we first try to find the image with the smallest attenuation
  // as given by the elixir keyword SKMAGATT (this keyword is not in all images)

  double min_attenuation = 12; 
  double attenuation;
  
  for (LightCurve::const_iterator it = begin(); it != end(); ++it) {
  //for (ReducedImageCIterator im=Images.begin(); im != Images.end(); ++im) {
    const Fiducial<PhotStar> *fs = *it;
#ifdef DEBUG
    cout << fs->Image()->FitsName() << " photomratio = " << fs->photomratio << endl;
#endif
    FitsHeader head(fs->Image()->FitsName());
    if(head.HasKey("SKMAGATT")) {
      attenuation = head.KeyVal("SKMAGATT");
      if(attenuation<min_attenuation) {
	min_attenuation = attenuation;
	photometric_image_fitsname = fs->Image()->FitsName();
	photomratio = fs->photomratio;
      }
    }
  }
  if(min_attenuation>10) {
    cout << "LightCurve::computeElixirZeroPoint WARNING no info on attenuation, using reference image zero point" << endl;
  }
  cout << "LightCurve::computeElixirZeroPoint min_attenuation= " << min_attenuation << endl;
  cout << "LightCurve::computeElixirZeroPoint photometric_image= " << photometric_image_fitsname << endl; 
  
#endif
  
  // we do not take into account color terms !!!
  FitsHeader refhead(photometric_image_fitsname);
  double expo =  refhead.KeyVal("TOADEXPO");
  double PHOT_C =  refhead.KeyVal("PHOT_C");
  double PHOT_K =  refhead.KeyVal("PHOT_K");
  double AIRMASS =  refhead.KeyVal("AIRMASS");
  
  return 2.5*log10(expo) + PHOT_C + PHOT_K*(AIRMASS-1.) + 2.5*log10(photomratio); // if photomratio>1, flux in ref < flux photometric => mag in ref > mag  photometric
}
