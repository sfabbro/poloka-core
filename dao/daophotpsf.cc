#include "cdaophot.h"
#include "daophot.h"
#include "daophotpsf.h"
#include "daophotio.h"


int DaoPsfVariability(const int nstars)
{
  if (nstars > 20) return 2;
  else if (nstars > 10) return 1;
  return 0;
}

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

void DaoPsf::Allocate(const int Npsf, const int Npar)
{
  if (table) delete [] table;
  if (param) delete [] param;
  table = new float[Npsf*Npsf*Npsf];
  param = new float[Npar];

}

bool DaoPsf::read(const string &FileName)
{
  if (!FileExists(FileName)) 
    {
      cerr << " DaoPsf::read() : psf file " << FileName << " not found \n";
      return false;
    }

  Allocate(MAXPSF, MAXPAR);

  int istat = RDPSF(FileName.c_str(), &type, param, &MAXPAR, &npar, 
		    table, &MAXPSF, &MAXEXP, &npsf, &nexp, &nfrac, 
		    &psfmag, &bright, &xpsf, &ypsf);

  if (istat == -1)
    {
      cerr << " DaoPsf::read() may have found cheese NaN " << endl;
      delete [] table;
      delete [] param;
      return false;
    }

  radius = (double(npsf-1)/2.-1.)/2.;

  return true;
}

double DaoPsf::ThetaXY() const
{
  if (type == 1) return 0;        // gaussian
  if (type < 5)  return param[2]; // lorentzians
  return param[3];                // pennys
}
		
void DaoPsf::dump(ostream &Stream) const
{
  size_t oldp = Stream.precision();
  Stream << " DaoPsf::dump() TYPE=" << Type()
	 << setiosflags(ios::fixed) << setprecision(3)
	 << " HWHMX=" << HwhmX() << " HWHMY=" << HwhmY()
	 << " THETAXY=" << ThetaXY() << endl;
  Stream << resetiosflags(ios::fixed) << setprecision(oldp);
  
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


double DaoPsf::Value(const int i, const int j, const double &Xc, const double &Yc, 
		     double &DpDx, double &DpDy) const
{
  // dx,dy are the distance from the center of the pixel to the center of the star
  float dx = i-Xc;
  float dy = j-Yc;

  // see a fortran example of use of usepsf in ../dao_stuff/addstar.f
  // there is no Xc-1 because Xc in in TOADS coordinate system (but xpsf is in daophot system)
  // deltax, deltay are relative frame coordinates
  float deltax = Xc/xpsf-1.;
  float deltay = Yc/ypsf-1.;

  // derivatives relatively to x and y
  float dvdxc, dvdyc;  

  double val = USEPSF(&type, &dx, &dy, &bright, param, table, 
		      &MAXPSF, &MAXPAR, &MAXEXP,
		      &npsf, &npar, &nexp, &nfrac, &deltax, &deltay, 
		      &dvdxc, &dvdyc);  

  // normalize daophot from internal cooking (see manual)
  double scale = pow(10, 0.4*(psfmag-DAOPHOT_APER_ZP));
  
  DpDx = dvdxc*scale;
  DpDy = dvdyc*scale;

  return val*scale;
}

/*
//! fit only one SEStar with a DaoPsf
void PeakFit(SEStar &Star, const DaoPsf& psf) const;
{
  float x      = Star.x - DAOPHOT_TOADS_SHIFT;
  float y      = Star.y - DAOPHOT_TOADS_SHIFT;
  float deltax = (x-1)/psf->xpsf -1 ;
  float deltay = (y-1)/psf->ypsf -1 ;
  float scale  = Star.flux;
  float sky    = Star.Fond();
  float perr   = opt[PercError].Value()*0.01;
  float pkerr  = opt[ProfError].Value()*0.01;
  float radius = min(opt[PsfRadius].Value(), float(((psf->npsf-1.)/2. - 1.)/2.) );
  int lx       = max(1, int(x-radius)+1);
  int ly       = max(1, int(y-radius)+1);
  int nx       = min(SIZE.ncol, int(x-radius)+1) - lx +1;
  int ny       = min(SIZE.nrow, int(y-radius)+1) - ly +1;
  x += -lx+1;
  y += -ly+1;
 
  const int maxbox = 69; // box of the litte vignette around the star
  float *f = new float[maxbox*maxbox];

  int status, niter;
  float errmag, chi, sharp;

  RDARAY("DATA", &lx, &ly, &nx, &ny, &maxbox, f, &status);

  PKFIT(f, &nx, &ny, &maxbox, &x, &y, &scale, &sky, &radius, 
	&lowbad, &opt[AduHighDatum].Value(), &opt[Gain].Value(), &opt[ReadNoise].Value(), &perr, &pkerr, 
	&psf->bright, &psf->type, psf->param, 
	&psf->MAXPAR, &psf->npar, psf->table, 
	&psf->MAXPSF, &psf->MAXEXP, &psf->npsf, &psf->nexp, 
	&psf->nfrac, &deltax, &deltay, &errmag, &chi, 
	&sharp, &niter, &global_sky);

  delete [] f;

  Star.x       = x + DAOPHOT_TOADS_SHIFT;
  Star.y       = y + DAOPHOT_TOADS_SHIFT;
  Star.flux    = scale;
  Star.EFlux() = errmag;
  Star.Chi()   = chi;
  Star.Sharp() = sharp;
  Star.Iter()  = niter;
}
*/
