#include "gtransfo.h"
#include "../dao_stuff/cdaophot.h"
#include "daophot.h"
#include "daophotpsf.h"
#include "daophotio.h"
#include "fileutils.h"
#include "imagematch.h"
#include "fastfinder.h"
#include "vutils.h"

//***************************** DaoPsf ****************************

DaoPsf::DaoPsf(const DbImage &DbIm) : param(0), table(0)
{
  if (!read(DbIm.ImagePsfName())) cerr << " DaoPsf::DaoPsf() not initialized" << endl;
}

DaoPsf::DaoPsf(const string &FileName) : param(0), table(0)
{
  if (!read(FileName)) cerr << " DaoPsf::DaoPsf() not initialized" << endl;
}

DaoPsf::~DaoPsf()
{
  if (table) delete [] table;
  if (param) delete [] param;
}

void DaoPsf::Allocate()
{
  if (table) delete [] table;
  if (param) delete [] param;
  table = new float[MAXPSF*MAXPSF*MAXPSF];
  param = new float[MAXPAR];

}

bool DaoPsf::read(const string &FileName)
{
  if (!FileExists(FileName)) 
    {
      cerr << " DaoPsf::read() : psf file " << FileName << " not found \n";
      return false;
    }

  Allocate();

  int istat = RDPSF(FileName.c_str(), &type, param, &MAXPAR, &npar, 
		    table, &MAXPSF, &MAXEXP, &npsf, &nexp, &nfrac, 
		    &psfmag, &bright, &xpsf, &ypsf);

  if (istat == -1)
    {
      cerr << " DaoPsf::read may have found cheese NaN " << endl;
      delete [] table;
      delete [] param;
      return false;
    }

  radius = (double(npsf-1)/2.-1.)/2.;

  return true;
}

double DaoPsf::ThetaXY() const
{
  if (type == 1) return 0;//gaussian
  if (type < 5) return param[2];//lorentzians
  return param[3];//pennys
}
		
void DaoPsf::dump(ostream &Stream) const
{
  Stream << " DaoPsf " << Type() << endl;
  Stream << " Parameters " << " HWHMX = " << HwhmX() << " HWHMY = " << HwhmY()
	 << " THETAXY = " << ThetaXY() << endl;
}

string DaoPsf::Type() const
{
  switch (type) 
    {
    case 1: return "GAUSSIAN";
    case 2: return "MOFFAT15";
    case 3: return "MOFFAT25";
    case 4: return "LORENTZ";
    case 5: return "PENNY1";
    case 6: return "PENNY2";
    }
  cerr << " DaoPsf::Type() unknown " << endl;
  return "UNKNOWN";
}

// returns value of psf on pixel (i,j) of the psf star centered on (Xc,Yc)
double DaoPsf::Value(const int i, const int j, const double &Xc, const double &Yc, 
		     double &DpDx, double &DpDy) const
{
  // dx,dy are the distance from the center of the pixel to the center of the star
  float dx = i-Xc;
  float dy = j-Yc;

  // see a fortran example of use of usepsf in dao_stuff/addstar.f
  // there is no Xc-1 because Xc in in TOADS coordinate system (but xpsf is in daophot system)
  // deltax, deltay are relative frame coordinates
  float deltax = Xc/xpsf-1;
  float deltay = Yc/ypsf-1;

  // derivatives relatively to x and y
  float dvdxc, dvdyc;
  
  double val = 0.;
  DpDx = 0;
  DpDy = 0;

  val = USEPSF(&type, &dx, &dy, &bright, param, table, 
	       &npsf, &npar, &nexp, &nfrac, &deltax, &deltay, 
	       &dvdxc, &dvdyc);  

  // normalize daophot from internal cooking (see manual)
  double scale = pow(10, 0.4*(psfmag-DAOPHOT_APER_ZP));
  
  DpDx = dvdxc*scale;
  DpDy = dvdyc*scale;
  return val*scale;
}



void DaoPsf::MakeStar(Kernel &ImStar, const BaseStar &Star) const
{
  if (ImStar.HSizeX() == 0 || ImStar.HSizeY() == 0) 
    ImStar = Kernel(int(radius), int(radius));
  int hx = ImStar.HSizeX();
  int hy = ImStar.HSizeY();
  int ic = int(Star.x);
  int jc = int(Star.y);
  for (int j=-hy; j<=hy; ++j) 
    for (int i=-hx; i<=hx; ++i)
      {
	double dpdx = 0 ,dpdy = 0;
	ImStar(i,j) = Star.flux * Value(ic+i,jc+j, Star.x, Star.y, dpdx, dpdy);
      }
}


bool UpdateSeeingFromDaoPsf(ReducedImage &Rim)
{
  DaoPsf psf;
  if (!psf.read(Rim.ImagePsfName()))
    {
      cerr << " UpdateFromPsf() : problem reading PSF file\n";
      return false;
    }
    
  cout << " Updating seeing from PSF parameters" << endl;
  psf.dump();
  // compute seeing a la delphine
  double sigmaX = psf.HwhmX() * 0.8493218;
  double sigmaY = psf.HwhmY() * 0.8493218;
  double thetaXY = psf.ThetaXY();
  // thetaXY is not from the gaussian.
  double r = (1-thetaXY*thetaXY);
  double seeing;
  if ((r > 1.e-10) && (sigmaX>1.e-10) && (sigmaY>1.e-10)) 
    seeing = sqrt(sqrt(r)*sigmaX*sigmaY);
  else seeing =  max(sigmaX,sigmaY); // is this ok?
  Rim.SetSeeing(seeing, "Seeing updated from DAOPHOT PSF file");
  Rim.SetPsfShapeParams(sigmaX, sigmaY, thetaXY, "Sigmas and theta updated from DAOPHOT PSF file");
  cout << " UpdateFromDaoPsf : Updated seeing = " << seeing << " sigmas" << endl;  
  return true;
}

//***************************** PsfConditions ****************************

PsfConditions::PsfConditions() :
  fw_min(0.5), fw_max(5.), 
  mx_min(1e-10), mx_max(15.), 
  my_min(1e-10), my_max(15.), 
  flux_min(1e-10), flux_max(1e30),
  fm_max(150000.), cstar_min(0.), 
  neigh_rad(10.), edges_rad(10.),
  flag_min(0), flag_max(0), flagbad_max(0)
{}

PsfConditions::PsfConditions(const ReducedImage &Rim) : 
  catalog(Rim.CatalogName()),
  mx_min(1e-10), mx_max(15.), 
  my_min(1e-10), my_max(15.),
  flux_min(1e-10), flux_max(1e38),  cstar_min(0.), 
  flag_min(0), flag_max(0), flagbad_max(0)
{
  double fwhm = Rim.Seeing()*2.3548;
  fw_min = fwhm * 0.3; fw_max = fwhm * 3;
  fm_max = Rim.Saturation();
  neigh_rad = fwhm*2; edges_rad = fwhm * 2;
  usableFrame = Rim.UsablePart();
}
  
bool PsfConditions::AreVerified(const SEStar* star) const 
{
  return (
	  (star->Fwhm() > fw_min && star->Fwhm() < fw_max) &&
	  (star->Mxx() > mx_min && star->Mxx() < mx_max)  &&
	  (star->Myy() > my_min && star->Myy() < my_max)  &&
	  (star->flux > flux_min && star->flux < flux_max) &&
	  (star->Flag() >= flag_min && star->Flag() <= flag_max) &&
	  (star->FlagBad() <= flagbad_max) &&
	  (star->Fluxmax() < fm_max) &&
	  (star->Cstar() > cstar_min) &&
	  (!catalog.HasCloseNeighbor(*star, neigh_rad)) &&
	  (usableFrame.MinDistToEdges(*star)> edges_rad)
	  );
}

//***************************** PsfStars ****************************


PsfStars::PsfStars(const ReducedImage &Rim) :SEStarList(Rim.CatalogName()), rim(&Rim)
{}

size_t PsfStars::FilterMultiIm(const ReducedImageList &RimList, const ReducedImage *ref, 
			    const double MaxDist) 
{
  for (ReducedImageCIterator rit=RimList.begin(); rit!=RimList.end(); ++rit)
    {
      const ReducedImage *rim = *rit;
      SEStarList currentCat(rim->CatalogName());
      Gtransfo *ref2cur=0, *cur2ref=0;
      ImageListMatch(*ref, *rim, ref2cur, cur2ref);
      FastFinder finder(*SE2Base(&currentCat));
      for (SEStarIterator it=begin(); it != end(); )
	{
	  Point refstar = ref2cur->apply(**it);
	  const BaseStar *neighbour = finder.FindClosest(refstar, MaxDist);
	  if (!neighbour) {++it;continue;}
	  double distance = refstar.Distance(*neighbour); 
	  if (distance < MaxDist) ++it;
	  else
	    {
	      it = erase(it);
	    }
	}
      delete ref2cur;
      delete cur2ref;
    }

  return size();
}

size_t PsfStars::FilterCond(const PsfConditions &Cond)
{
  // get rid of surely bad objects from the list
  for (SEStarIterator it=begin(); it != end(); )
    {      
      if (Cond.AreVerified(*it)) ++it;
      else
	{
	  it = erase(it);
	}  	 	
    }

  // sort per brightness
  sort(&DecFluxMax);
  
  return size();				   
}

bool DecreasingChi(const SEStar *S1, const SEStar *S2)
{
  return (S1->Chi() > S2->Chi());
}

size_t PsfStars::FilterNei()
{

  // read neighbor file
  SEStarList neiList;
  string neifile = rim->ImageCatalogName(DaophotNei);
  read_dao<DaophotNei>(neifile, neiList);
  neiList.sort(&DecreasingChi);
  SEStarCIterator itn = neiList.begin(); 
  double cutchi;
  while (itn != neiList.end() && (*itn)->Chi() > 0) ++itn;
  cutchi = 2*(*--itn)->Chi();
  cout << " cutchi = " << cutchi << endl;

  cout << " filter nei before " << size() << endl;
  
  // loop on PSF stars. 
  SEStarCIterator endnei = neiList.end();
  for (SEStarIterator it=begin(); it != end(); )
    {
      // loop on neighbor stars
      SEStarCIterator itnei = neiList.begin();
      while ((itnei != endnei)&& ((*itnei)->N() != (*it)->N())) ++itnei;

      // cut on unfound stars, flagged, and bad chis
      if ((itnei != endnei) && 
	  ((*itnei)->Flag() < 256)  &&
	  ((*itnei)->Chi() > 0.) &&
	  ((*itnei)->Chi() < cutchi) &&
	  ((*itnei)->Chi() < 0.1)) 
	++it;
      else
	{
	  it = erase(it);
	}
    }

  cout << " filter nei after  " << size() << endl;

  return size();
}


size_t PsfStars::FilterAls(const double ChiMax, 
			   const double SharpMin, 
			   const double SharpMax)
{
  // read the allstar file
  SEStarList alsList;
  string alsfile = rim->ImageCatalogName(DaophotAls);
  read_dao<DaophotAls>(alsfile, alsList);
  
  // cout << " filter als before " << size() << endl;

  // loop on PSF stars. 
  SEStarIterator endals = alsList.end();
  for (SEStarIterator it = begin(); it != end(); )
    {
      // loop on als stars. 
      SEStarCIterator itals = alsList.begin();
      while ((itals != endals) && ((*itals)->N() != (*it)->N())) ++itals;
      // follow Stetson instructions (DAOPHOT II manual): 
      // cuts on sharp and chi derived from ALLSTAR
      if ((itals != endals)  &&
	  ((*itals)->Chi() < ChiMax) &&
	  ((*itals)->Sharp() < SharpMax) &&
	  ((*itals)->Sharp() > SharpMin))
	++it;
      else
	{
	  it = erase(it);
	}
    }   

  // cout << " filter als after  " << size() << endl;

  return size();
}

