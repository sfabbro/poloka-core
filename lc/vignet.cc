#include <iomanip>
#include "vignet.h"

//#define DEBUG

void Vignet::Allocate(const int Nx, const int Ny)
{
  Data.Allocate(Nx,Ny);
  Weight.Allocate(Nx,Ny);
  Resid.Allocate(Nx,Ny);
}

bool Vignet::Load(const PhotStar *AStar)
{    
#ifdef DEBUG
  cout << " Vignet::Load(" << AStar<< ");" << endl;
  cout << " rim->refCount() = " << rim->refCount() << endl;
  cout << " rim->Name() = " << rim->Name() << endl;
#endif

  if (!AStar) return false;

  Star = AStar;

  int xc = int(Star->x);
  int yc = int(Star->y);
  
  xstart = max(0, xc-Data.HSizeX());
  ystart = max(0, yc-Data.HSizeY());
  xend   = min(xc+Data.HSizeX()+1, rim->XSize());
  yend   = min(yc+Data.HSizeY()+1, rim->YSize());
  
  if (!rim->HasImage())
    {
      cerr << " Vignet::Load() : Error : " << rim->Name() << " does not have image" << endl;
      return false;
    } 
  
  Data.readFromImage(rim->FitsName(), *this);

  if (rim->HasWeight()) Weight.readFromImage(rim->FitsWeightName(), *this);
  
#ifdef DEBUG
  cout << " in Vignet::Load Data.HSizeX() = " << Data.HSizeX() << endl;
  cout << " in Vignet::Load Data.HSizeY() = " << Data.HSizeY() << endl;
#endif
  return true;
}


void Vignet::Resize(const int HNewX, const int HNewY)
{
#ifdef DEBUG
  cout << " Vignet::Resize(" <<  HNewX << "," << HNewY << "); rim is " << rim->Name() << endl;
#endif

  if  (HNewX < 0 || HNewY < 0) 
    {
      cerr << " Vignet::Resize(" << HNewX << "," << HNewY << ") : impossible \n";
      return;
    } 

  hx = HNewX;
  hy = HNewY;

  Allocate(2*hx+1,2*hy+1);
#ifdef DEBUG
  cout << " in Vignet::Resize Data.HSizeX() = " << Data.HSizeX() << endl;
  cout << " in Vignet::Resize Data.HSizeY() = " << Data.HSizeY() << endl;
#endif
}

void Vignet::Resize(const double &ScaleFactor)
{
#ifdef DEBUG
  cout << " Vignet::Resize(" << ScaleFactor << ");" << endl;
#endif
  
  if (ScaleFactor*hx > Data.HSizeX() || ScaleFactor*hy > Data.HSizeY() || ScaleFactor < 0.)
    {
      cerr << " Vignet::Resize(" << ScaleFactor << ") : Error:  can't resize, max scale is " 
	   << max(double(Data.HSizeX())/hx, double(Data.HSizeY())/hy) << endl;
      return;
    }
  
  hx = int(ceil(ScaleFactor*double(hx)));
  hy = int(ceil(ScaleFactor*double(hy)));

  Allocate(2*hx+1,2*hy+1);
}

bool Vignet::IsInside(const Point& point) {
  cerr << "Vignet::IsInside TO MODIFY !!!!!!!!!!!!!!" << endl;
  // JG pas sur de lui
  return (point.x >= xstart) 
    && (point.x < xend)
    && (point.y >= ystart)
    && (point.y < yend);
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

