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
  
  // we need this check otherwise we get overflow in USEPSF
  if(fabs(dx)>radius || fabs(dy)>radius) {
    //cout << "return 0" << endl;
    DpDx=0;
    DpDy=0;
    return 0;
  }

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

  //if(!((val*scale)>0)) {
  // cout << "DaoPsf::Value val,psfmag,scale " << val << "," << psfmag << "," << scale << endl;
  //}

  return val*scale;
}

double DaoPsf::Value(const int i, const int j, const Point &Pt,
		     double &DpDx, double &DpDy) const
{
  return Value(i,j,Pt.x,Pt.y,DpDx,DpDy);
}

double DaoPsf::Value(const int i, const int j, const double &Xc, const double &Yc) const
{
  // derivatives relatively to x and y
  double dpdx, dpdy;  
  return Value(i,j,Xc,Yc,dpdx,dpdy);
}

double DaoPsf::Value(const int i, const int j, const Point &Pt) const
{
  return Value(i,j,Pt.x,Pt.y);
}
