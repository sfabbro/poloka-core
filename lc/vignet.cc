#include <iomanip>
#include "vignet.h"

bool Vignet::Load(const PhotStar *AStar)
{    
  if (!AStar) return false;

  cout << " loading " << *AStar << endl;

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
  
  return true;
}

void Vignet::Resize(const int HNewX, const int HNewY)
{
  if ((HNewX > Data.HSizeX()) ||(HNewY > Data.HSizeY()) || (HNewX < 0) ||(HNewY < 0))
    {
      cerr << " Vignet::Resize(" << HNewX << "," << HNewY << ") : Error :  can't resize, max half sizes are (" 
	   << Data.HSizeX() << "," << Data.HSizeY() << ")\n";
      return;
    }
  
  hx = HNewX;
  hy = HNewY;
  
}

void Vignet::Resize(const double &ScaleFactor)
{
  
  if (ScaleFactor*hx > Data.HSizeX() || ScaleFactor*hy > Data.HSizeY() || ScaleFactor < 0.)
    {
      cerr << " Vignet::Resize(" << ScaleFactor << ") : Error:  can't resize, max scale is " 
	   << max(double(Data.HSizeX())/hx, double(Data.HSizeY())/hy) << endl;
      return;
    }
  
  hx = int(ceil(ScaleFactor*double(hx)));
  hy = int(ceil(ScaleFactor*double(hy)));
}

bool Vignet::IsInside(const Point& Pt) const
{
  return (Pt.x > xstart) && (Pt.x < xend) 
    &&   (Pt.y > ystart) && (Pt.y < yend);
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
    
  stream << " Star " << *Vig.Star << endl;
  
  return stream;
}

