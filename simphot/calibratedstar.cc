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
  u = 99 ;
  if ( D.HasKey("mu") )
    u = D.getval("mu");
  g = D.getval("mg");
  r = D.getval("mr");
  i = D.getval("mi");  
  z = D.getval("mz"); 
  ue = 99 ;
  if ( D.HasKey("emu") ) 
    ue = D.getval("emu");
  ge = D.getval("emg");
  re = D.getval("emr");
  ie = D.getval("emi");
  ze = D.getval("emz");
  flux = D.getval("flux"); // what can it be??
  id = D.Rank();
  neighborDist = -1 ; neighborFlux = -1 ;
  neighborFluxContamination = -1 ;
  neighborNsigma = -1 ;
}


#include "imageutils.h"
static Frame FrameApplyTransfo(const Frame& InputFrame, const Gtransfo &T, 
			       const WhichTransformed W)
{
  return ApplyTransfo(InputFrame, T, W);
}

#include "gtransfo.h"
#include "frame.h"
#include "wcsutils.h"
#include "fitsimage.h"

#include "polokaexception.h"
#include "sstream"

#include "fastfinder.h"
#include "apersestar.h"

static CalibratedStar * Checked_CalibratedStar(FastFinder const & finder, const Point & pix, const DicStar & ds_modele, double reference_seeing )
{
  const BaseStar * closest_basestar = NULL;  
  const AperSEStar * second_closest_star 
    = dynamic_cast<const AperSEStar *>(finder.SecondClosest(pix,50.,closest_basestar)); 
  if ( !closest_basestar) {
    cerr << "ERROR cannot find star from photometric catalog at " << pix.x << " " << pix.y << endl;
    return NULL ;
  }

  CalibratedStar *cs = new CalibratedStar(ds_modele); // dc recopiee : cs a dc le meme id
  // mod coordinates
  cs->x = closest_basestar->x;
  cs->y = closest_basestar->y;

  // compute rough flux contamination
  const AperSEStar * closest_star = dynamic_cast<const AperSEStar *>(closest_basestar);
  double dist = 0;
  double max_flux_contamination = 0;
  double nsigma = 1.e6;

  // based on D2 ccd_35 : average is 2.378, 90%<3.1
  double calibration_seeing = 3.1;
  double calibration_aperture_radius = 5.*calibration_seeing;
  double seeing_scale = (calibration_seeing/reference_seeing);
  seeing_scale *= seeing_scale ;
  if ( second_closest_star ) {
      
    const double &aper_radius = calibration_aperture_radius;
    double mxx = second_closest_star->gmxx;
    double myy = second_closest_star->gmyy;
    double mxy = second_closest_star->gmxy;

    // scale with best image seeing
    mxx *= seeing_scale;
    myy *= seeing_scale;
    mxy *= seeing_scale;
      
      
    Point dP =  *second_closest_star - *closest_basestar;
    dist = sqrt(dP.x*dP.x+dP.y*dP.y);
      
    nsigma = 0;
    bool is_inside = (dist<aper_radius);
    double distance_scale = (dist - aper_radius)/dist;
    if(dist<aper_radius) {
      //cout << "whao , second_nearest position is inside aperture radius" << endl;
      distance_scale *= -1;
    }
	
    dP.x *= distance_scale;
    dP.y *= distance_scale;
      
    // calcul du moment de l'objet   
    double det = mxx*myy-mxy*mxy;
    nsigma = sqrt((myy*dP.x*dP.x-2.*mxy*dP.x*dP.y+mxx*dP.y*dP.y)/det);
    // for gaussian
    double frac_gaus= erfc(nsigma/sqrt(2.))/2.;
      
    // for a moffat 2.5
    // ------------------------------------------
    double nsigma_moffat = nsigma*0.6547; // inside ~4 sigma radius 
    double nsigma2 = nsigma_moffat*nsigma_moffat;
    double frac_moffat = 0.1591549431*(-2*nsigma_moffat+(3.141592654-2*atan(nsigma_moffat))*(1.+nsigma2))/(1.+nsigma2);
      
    if(is_inside) {
      frac_gaus   = 1.-frac_gaus;
      frac_moffat = 1.-frac_moffat;
    }

    if(frac_moffat>frac_gaus)
      max_flux_contamination = frac_moffat*second_closest_star->flux;
    else
      max_flux_contamination = frac_gaus*second_closest_star->flux;
      
    cout << "x,y,dist,ndist,mxx,mxy,myy,nsigma,fg,fm,rfrac="
	 << closest_star->x << " "
	 << closest_star->y << " " 
	 << dist << " " 
	 << dist - aper_radius << " " 	
	 << mxx  << " " 
	 << mxy  << " " 
	 << myy  << " " 
	 << nsigma  << " " 
	 << frac_gaus << " " 
	 << frac_moffat << " " 
	 << max_flux_contamination/closest_star->flux << " "
	 << endl;
  } // if second_closest
    
  if(second_closest_star) {
    cs->neighborDist = dist;
    cs->neighborFlux = second_closest_star->flux;
      }else{
	cs->neighborDist = closest_star->neighborDist;
	cs->neighborFlux = closest_star->neighborFlux;
      }
  cs->neighborFluxContamination = max_flux_contamination;
  cs->neighborNsigma = nsigma;

  return(cs);
}




static string getkey(const string& prefix, const DicStarList& catalog) {
  if(!catalog.HasKey(prefix)) {
    cerr << "catalog does not not have info about " << prefix << endl;
    exit(2);
    return "absent";
  }
  return prefix;
}



CalibratedStarList::CalibratedStarList(const string &CatalogName, 
				       const  ReducedImageRef & refimage)
{


  FitsHeader geomHead(refimage->FitsName());
  GtransfoRef wcs = WCSFromHeader(geomHead);
  if (!wcs)
    {
      cout << "LightCurveFile::SimPhotFitForCalib : could not get wcs from image " << geomHead.FileName() << endl;
      exit(-1);
    }
  Frame imageFrame(geomHead);
  imageFrame.CutMargin(-100); //actually enlarges the frame


  DicStarList dsl(CatalogName);
  // check that the image and catalog more or less match
  Frame radecImageFrame = FrameApplyTransfo(imageFrame, *wcs, LargeFrame);
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

  Gtransfo *radec2pix = wcs->InverseTransfo(0.01, imageFrame);

  // band
  double mag_limit = 99 ;
  string band = refimage->Band();
  string mag_key=getkey("m"+band,dsl);
  
  if(band=="g") mag_limit = 21.;
  if(band=="r") mag_limit = 21.;
  if(band=="i") mag_limit = 21.;
  if(band=="z") mag_limit = 21.;
  


  // apersestarcat for ref
  string refimage_catalog = refimage->AperCatalogName();
  AperSEStarList starlist(refimage_catalog);
  FastFinder finder((const BaseStarList&)starlist);

  double mag ;
  for (DicStarCIterator i = dsl.begin(); i != dsl.end(); ++i)
    {
      const DicStar &ds =  **i;

      mag=ds.getval(mag_key);
      if(mag>mag_limit) continue; // ignore dim stars unused for calibration to save CPU


      Point radec(ds);
      Point pix = radec2pix->apply(radec);
      if (!imageFrame.InFrame(pix)) continue;

      // check for closest and second closest in aperselist to : get better coordinates and compute contamination
      double reference_seeing = refimage->Seeing();
      CalibratedStar *cs = Checked_CalibratedStar(finder,pix, ds,reference_seeing );
      if (cs ) 
	{
	  cs->flux = pow(10.,-0.4*mag);
	  push_back(cs); 
	} 
    }
  delete radec2pix;
}

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

