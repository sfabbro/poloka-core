#include "calibratedstar.h"

#include "dicstar.h"


/* tags of the catalogs :
# x:
# y:
# flux:
# mu:
# emu:
# nmesu:
# dccdu:
# mg:
# emg:
# nmesg:
# dccdg:
# mr:
# emr:
# nmesr:
# dccdr:
# mi:
# emi:
# nmesi:
# dccdi:
# mz:
# emz:
# nmesz:
# dccdz:
# level:
# mu_ab:
# mg_ab:
# mr_ab:
# mi_ab:
# mz_ab:
# mu_ab:
# mg_ab:
# mr_ab:
# mi_ab:
# mz_ab:
# end

*/


CalibratedStar::CalibratedStar(const DicStar &D) : BaseStar(D)
{
  ra = D.x;
  dec = D.y;
  u = D.getval("mu");
  g = D.getval("mg");
  r = D.getval("mr");
  i = D.getval("mi");  
  z = D.getval("mz");  
  ue = D.getval("emu");
  ge = D.getval("emg");
  re = D.getval("emr");
  ie = D.getval("emi");
  ze = D.getval("emz");
  flux = D.getval("flux"); // what can it be??
  id = D.Rank();
}


#include "imageutils.h"
static Frame FrameApplyTransfo(const Frame& InputFrame, const Gtransfo &T, 
			       const WhichTransformed W)
{
  return ApplyTransfo(InputFrame, T, W);
}

#include "gtransfo.h"
#include "frame.h"

#include "polokaexception.h"
#include "sstream"


CalibratedStarList::CalibratedStarList(const string &CatalogName, const Gtransfo *WCS, const Frame& ImageFrame)
{
  DicStarList dsl(CatalogName);
  // check that the image and catalog more or less match
  Frame radecImageFrame = FrameApplyTransfo(ImageFrame, *WCS, LargeFrame);
  DicStarList temp;
  dsl.ExtractInFrame(temp, radecImageFrame);
  if (temp.size()==0)
    {
      stringstream mess;
      mess << " don't not find any overlap between the catalog " 
	   << CatalogName << " and the geom ref image " 
	     << " we stop here " << endl;
      cout << mess;
      throw (PolokaException(mess.str()));
    }

  Gtransfo *radec2pix = WCS->InverseTransfo(0.01, ImageFrame);
  for (DicStarCIterator i = dsl.begin(); i != dsl.end(); ++i)
    {
      const DicStar &ds =  **i;
      Point radec(ds);
      Point pix = radec2pix->apply(radec);
      if (!ImageFrame.InFrame(pix)) continue;
      CalibratedStar *cs = new CalibratedStar(ds);
      cs->x = pix.x;
      cs->y = pix.y;
      push_back(cs);      
    }
  delete radec2pix;
}

