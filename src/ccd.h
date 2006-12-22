#ifndef CCD__H 
#define CCD__H



#include "fitsimage.h"
#include "wcsutils.h"
#include "reducedimage.h"
#include "gtransfo.h"
#include "imageutils.h"
namespace Toads_ccd{
class Ccd {
  private :

  ReducedImageRef ri;
  Frame whereInRef;
  bool keepIt;

public :

  Ccd(const std::string &ReducedImageName, const Gtransfo &RefRaDec2Pix) 
    : ri(new ReducedImage(ReducedImageName))
  {
    keepIt = false;
    if (!ri->IsValid()) 
      {
	std::cerr << " cannot find " << ReducedImageName << std::endl;
	ri = ReducedImageRef((ReducedImage *)NULL);
	return;
      }
    Gtransfo *wcs;
    FitsHeader head(ri->FitsName());
    WCSFromHeader(head,wcs);
    Gtransfo *pix2RaDec = dynamic_cast<TanPix2RaDec*>(wcs);
    if (!pix2RaDec)
      {
	std::cerr << " no tan wcs in image " << Name() 
		  << " : we ignore it" << std::endl;
	return;
      }
    if (!head.HasKey("ZP"))
      {
	std::cerr << " no ZP key in " << Name() 
		  << " : we ignore it" << std::endl;
	return;
      }	
    Gtransfo* pix2RefPix  = GtransfoCompose(&RefRaDec2Pix, pix2RaDec);
    whereInRef = ApplyTransfo(ri->UsablePart(),*pix2RefPix, LargeFrame);
    delete wcs;
    keepIt = true;
  }

  std::string Name() const { return ri->Name();}

  bool IsValid() const { return keepIt && ri && ri->IsValid();}

  Frame WhereInRef() const { return whereInRef;}

};

class Ccds : vector<Ccd>
{
private :
  Gtransfo *raDec2RefPix;
  std::string refName;
  Frame wholeRefFrame;

public :
  Ccds(const ReducedImage &Ref)
  {
    raDec2RefPix = NULL;
    refName = Ref.Name();
    FitsHeader head(Ref.FitsName());
    Gtransfo *wcs;
    if (!WCSFromHeader(head,wcs))
      {
	std::cerr << " no WCS in " << refName << " giving up " << std::endl;
	return;
      }
    wholeRefFrame = Frame(head);
    raDec2RefPix = wcs->InverseTransfo(0.1 /*precision in pixels*/, 
				       wholeRefFrame);
    delete wcs;
  }

  const Frame &WholeRefFrame() const { return wholeRefFrame;}

  const Ccd *LocateImage(const std::string &Name)
  {
    for (unsigned k=0; k < size(); ++k)
      {
	Ccd *ccd = &(*this)[k];
	if (ccd->Name() == Name) return ccd;
      }
    return NULL;
  }

  bool AddImage(const std::string &Name)
  {
    Ccd ccd(Name, *raDec2RefPix);
    if (ccd.IsValid()) { push_back(ccd); return true;}
    else return false;
  }

  void SelectInsiders(const Frame &SubFrame, StringList &Insiders)
  {
    Insiders.clear();
    for (unsigned k=0; k < size(); ++k)
      {
	Ccd &ccd = (*this)[k];
	Frame intersection = SubFrame*ccd.WhereInRef();
	if (intersection.Area() >0) Insiders.push_back(ccd.Name());
      }
  }
    

};};
#endif
