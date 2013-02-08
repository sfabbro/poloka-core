#include <iostream>
#include <cstring> /* for memset */

#include <poloka/histo1d.h>

using namespace std;

Histo1d::Histo1d(int Nx, float Minx, float Maxx) : nx(Nx), minx(Minx)
{
if (Maxx != minx)
  scalex = Nx/(Maxx - minx);
else 
  {
  cerr << " Histo1d :: minx = maxx requested" << endl;
  scalex = 1;
  }
data = new float[nx];
memset(data,0,nx*sizeof(float));
}

double Histo1d::MaxBin(double &X) const
{
  float *p ;
  int imax=0;
  float valmax = -1e30;
  p = data ;
  for (int ibin = 0 ; ibin < nx ; ibin++, p++ )
    {
      if (*p > valmax) {valmax = *p; imax = ibin;}
    }
  //X = minx + ((float)imax + 0.5)/scalex;
  X = minx + ((float)imax )/scalex;
  return valmax;
}

void Histo1d::ZeroBins(double Xmin, double Xmax)
{
  int binmin = int(floor((Xmin - minx)*scalex));
  int binmax = int(floor((Xmax - minx)*scalex));
  for (int bin = binmin; bin < binmax; bin ++)
    {
      if (bin<0 || bin>=nx) return; 
      data[bin] = 0;
    }
}

int Histo1d::BinAt(double X) const { 
  int bin = int(floor((X - minx)*scalex));
  if (bin<0 || bin>=nx) return -1;
  return bin;
} 

double Histo1d::BinCenter(int bin) const{
  if (bin<0 || bin>=nx) {
    cerr << "Histo1d::BinCenter ERROR out of range" << endl;
    return 0;
  }
  return (minx + (bin+0.5)/scalex);
}

void Histo1d::dump(ostream &stream) const {
  stream << nx << " " << minx << " " << minx+nx/scalex << endl;
  for(int bin=0;bin<nx;++bin)
    stream << bin << " " << BinCenter(bin) << " " << data[bin] << endl;
}
