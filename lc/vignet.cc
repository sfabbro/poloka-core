#include <iomanip>
#include "vignet.h"

#define FNAME
#define DEBUG

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


void Vignet::Resize(const int Hx, const int Hy)
{
#ifdef FNAME
  cout << " > Vignet::Resize(const int Hx, const int Hy) Hx,Hy = " << Hx << "," << Hy << endl;
#endif

  if  (Hx < 0 || Hy < 0) 
    {
      cerr << " Vignet::Resize(" << Hx << "," << Hy << ") : Error : impossible \n";
      return;
    } 
  if (!Star) {
    cerr << " Vignet::Resize : Error no star loaded " << endl;
  }
  
  hx=Hx;
  hy=Hy;
  
  Load(Star);
  
#ifdef DEBUG
  cout << " in Vignet::Resize (2) Data.HSizeX(),Data.HSizeY() = " << Data.HSizeX() << "," << Data.HSizeY()<< endl;
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
  DPixel *pw, *pres;
  for (int j=-hy; j<=hy; ++j)
    {
      pw   = &Weight(-hx,j);
      pres = &Resid (-hx,j);
	for (int i=-hx; i<=hx; ++i, ++pres) 
	  chi2 += *pw++ * (*pres) * (*pres) ;
    }
  
  return chi2;
}

double Vignet::MeanResid() const
{
  
  double mean = 0.;
  double npix = 0.;
  DPixel *pres;
  for (int j=-hy; j<=hy; ++j)
    {
      pres = &Resid(-hx,j);
	for (int i=-hx; i<=hx; ++i) 
	  {
	    mean += (*pres++);
	    npix++;
	  }
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
  
  if (!Vig.Resid.IsEmpty())
    stream << " chi2 = " << Vig.Chi2() 
	   << " mean resid = " << Vig.MeanResid() << endl;
    
  if (Vig.Star) stream << " Star " << *Vig.Star << endl;
  
  return stream;
}

