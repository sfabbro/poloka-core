#include <iostream>
#include "histo1d.h"
#include <string.h> /* for memset */

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

double Histo1d::MaxBin(double &X)
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
