#include <iomanip>
#include "vignet.h"

void Vignet::Allocate()
{
#ifdef FNAME
  cout << " > Vignet::Allocate() Nx,Ny = " << Nx() << "," << Ny() << endl;
#endif
  Data.Allocate(Nx(),Ny());
  Weight.Allocate(Nx(),Ny());
  Resid.Allocate(Nx(),Ny());
}

bool Vignet::Load(const PhotStar *AStar)
{    
#ifdef FNAME
  cout << " > Vignet::Load(const PhotStar *AStar) AStar=" << AStar << endl;
#endif

#ifdef DEBUG
  cout << "   rim->refCount() = " << rim->refCount() << endl;
  cout << "   rim->Name() = " << rim->Name() << endl;
#endif

  if (!AStar) return false;

  Star = AStar;

  int xc = int(Star->x);
  int yc = int(Star->y);
  
  xstart = max(0, xc-hx);
  ystart = max(0, yc-hy);
  xend   = min(xc+hx+1, rim->XSize());
  yend   = min(yc+hy+1, rim->YSize());
  
  Allocate();

  if (!rim->HasImage())
    {
      cerr << " Vignet::Load() : Error : " << rim->Name() << " does not have image" << endl;
      return false;
    } 
  
  Data.readFromImage(rim->FitsName(), *this);

  if (rim->HasWeight()) Weight.readFromImage(rim->FitsWeightName(), *this);
  
  return true;
}


void Vignet::Resize(const int Hx_new, const int Hy_new)
{
#ifdef FNAME
  cout << " > Vignet::Resize(const int Hx_new, const int Hy_new) Hx_new,Hy_new = " << Hx_new << "," << Hy_new << endl;
#endif

  if  (Hx_new < 0 || Hy_new < 0) 
    {
      cerr << " Vignet::Resize(" << Hx_new << "," << Hy_new << ") : Error : impossible \n";
      return;
    } 
  if (!Star) {
    cerr << " Vignet::Resize : Error no star loaded " << endl;
  }
  
  hx=Hx_new;
  hy=Hy_new;
  
  Load(Star);
  
#ifdef DEBUG
  cout << "  in Vignet::Resize, Data.HSizeX(),Data.HSizeY() = " << Data.HSizeX() << "," << Data.HSizeY()<< endl;
  cout << "  in Vignet::Resize, Hx(),Hy() = " << Hx() << "," << Hy()<< endl;
#endif
  
}

void Vignet::Resize(const double &ScaleFactor)
{

#ifdef DEBUG
  cout << " Vignet::Resize(" << ScaleFactor << ");" << endl;
#endif
 
 if (ScaleFactor <0)
    {
      cerr << " Vignet::Resize(" << ScaleFactor << ") : Error: impossible \n";
      return;
    }
  
 Resize(int(ceil(ScaleFactor*double(hx))),
	int(ceil(ScaleFactor*double(hy))));
}

bool Vignet::ShiftCenter(const Point& Shift)
{
  // check if shift does not go a chaille
  if (!IsInside(*Star+Shift))
    {
      cerr << " Vignet::ShiftCenter() : Error : " << Shift << " shift Star center outside of vignet \n";
      return false;
    }
  
  cout << "Vignet::ShiftCenter() : " << Shift  << endl;
  
  if (Star)
    {
      Star->x += Shift.x;
      Star->y += Shift.y;
    }

  return true;
}

void Vignet::ClearResidZeroWeight() {
  if (Resid.IsEmpty() || Weight.IsEmpty()) return;
  DPixel *pw, *pres;
  pw=Weight.begin();
  pres=Resid.begin();
  for (int i=0; i< Nx()*Ny(); ++i , ++pres, ++pw) {
    if(*pw<1.e-30)
      *pres=0;
  }
}

double Vignet::SigmaResid() const {
  if (Resid.IsEmpty() || Weight.IsEmpty()) return 0;
  DPixel *pw, *pres;
  pw=Weight.begin();
  pres=Resid.begin();
  double sw = 0;
  double sf = 0;
  double sf2 = 0;
  for (int i=0; i< Nx()*Ny(); ++i , ++pres, ++pw) {
    sw += *pw;
    sf += *pw * (*pres);
    sf2 += *pw * (*pres) * (*pres);
  }
  if(sw<1.e-30)
    return 0;
  return sqrt(sf2/sw-pow(sf/sw,2));
}

int  Vignet::NValidPixels() const {
  if (Resid.IsEmpty() || Weight.IsEmpty()) return 0;
  DPixel *pw = Weight.begin();
  int nok = 0;
  for (int i=0; i< Nx()*Ny(); ++i , ++pw) {
    if(*pw>1e-30)
      nok++;
  }
  return nok;
}


double Vignet::MaxPixResid() const {
  if (Resid.IsEmpty() || Weight.IsEmpty()) return 0;
  DPixel *pw, *pres;
  pw=Weight.begin();
  pres=Resid.begin();
  double max = 0;
  for (int i=0; i< Nx()*Ny(); ++i , ++pres, ++pw) {
    if(*pw>1.e-30) {
      if(fabs(*pres)>max)
	max = fabs(*pres);
    }
  }
  return max;
}



void Vignet::KillOutliers(const double& nsigma) {
 
  if (Resid.IsEmpty() || Weight.IsEmpty()) return;
  cout << " > Vignet::KillOutliers with nsigma = " << nsigma << endl;
  DPixel *pw, *pres;
  
  double sw = 0;
  double sf = 0;
  double sf2 = 0;
  pw=Weight.begin();
  pres=Resid.begin();
  for (int i=0; i< Nx()*Ny(); ++i , ++pres, ++pw) {
    sw += *pw;
    sf += *pw * (*pres);
    sf2 += *pw * (*pres) * (*pres);
  }
  double mean = sf/sw;
  double sigma = sqrt(sf2/sw-mean*mean);
  double thres = nsigma*sigma;
  pw=Weight.begin();
  pres=Resid.begin();
  int nbad = 0;
  for (int i=0; i< Nx()*Ny(); ++i , ++pres, ++pw) {
    if (fabs(*pres-mean)>thres) {
      *pw = 0;
      nbad ++;
    }
  }
  cout << "   in Vignet::KillOutliers, nbad,mean,sigma = " << nbad << ","  << mean << ","  << sigma << endl;
}

void Vignet::RobustifyWeight(const double& alpha, const double& beta)
{
  if (Resid.IsEmpty() || Weight.IsEmpty()) return;
  
  cout << " Vignet::RobustifyWeight(" << alpha << "," << beta << ") : re-weight" << endl;
  
  DPixel *pw, *pres;

  for (int j=-hy; j<=hy; ++j)
    {
      pw   = &Weight(-hx,j);
      pres = &Resid (-hx,j);
      for (int i=-hx; i<=hx; ++i, ++pw, ++pres)
	if (*pw != 0) 
	  *pw *=  1. / (1. + pow(fabs(*pres)/(sqrt(1./ *pw)*alpha), beta));
    }
}
  
double Vignet::Chi2() const
{
  if (Resid.IsEmpty() || Weight.IsEmpty()) return -1.;

  double chi2 = 0.;
  DPixel *pw=Weight.begin(), *pres=Resid.begin();
  for (int i=Nx()*Ny(); i; --i, ++pres)
    chi2 += *pw++ * (*pres) * (*pres);
  
  return chi2;
}

double Vignet::MeanResid() const
{
  
  double mean = 0.;
  double npix = 0.;
  DPixel *pres = Resid.begin();
  for (int i=Nx()*Ny(); i; --i)
    {
      mean += (*pres++);
      npix++;
    }
  if (npix != 0.) return mean/npix;

  return 0.;
}

  
void Vignet::writeInImage(const string& FitsFileName) const 
{ 
  Data.writeInImage(FitsFileName, *this); 
}


void writeModel(const Vignet& Vig, const string& FitsFileName)
{  
  Kernel model(Vig.Data);
  model -= Vig.Resid;
  model.writeInImage(FitsFileName, Vig);
}
  
ostream& operator << (ostream & stream, const Vignet& Vig)
{
  stream << " hx = " << Vig.hx 
	 << " hy = " << Vig.hy
	 << " rect = " << (Window) Vig << endl;
  
  if (!Vig.Resid.IsEmpty()) {
    stream << " Resid chi2 = " << Vig.Chi2() 
	   << " mean resid = " << Vig.MeanResid()
	   << " sigma resid = " << Vig.SigmaResid()
	   << " max pix resid = " << Vig.MaxPixResid()
	   << " valid pixels = " << Vig.NValidPixels()
	   << endl;
  }
  if (Vig.rim)  stream << " Image " << Vig.rim->Name() << endl;
  if (Vig.Star) stream << " Star " << *Vig.Star << endl;
  
  return stream;
}

