#include <iostream>
#include "histo2d.h"
#include <math.h> /* for floor */
#include <string.h> /* for memset*/

using namespace std;

Histo2d::Histo2d(int nnx, float mminx, float mmaxx, int nny, float mminy,float mmaxy)
{
nx = nnx;
ny = nny;
minx = mminx;
miny = mminy;
if (mmaxx!= mminx) 
  scalex = nx/(mmaxx-mminx); 
else 
  {
  cerr << " Histo2d: minx = maxx requested" << endl;
  scalex = 1.0;
  }
if (mmaxy != mminy)
  scaley = ny/(mmaxy-mminy);
else
  {
  cerr << " Histo2d : maxy = miny requested" << endl;
  scaley = 1.0;
  }
data = new float[nx*ny];
memset(data, 0, nx*ny*sizeof(float));
}

void Histo2d::Fill(float X, float Y, float Weight)
{
int ix, iy;
ix = (int) floor(( X - minx)*scalex);
if (ix <0 || ix >= nx) return;
iy = (int) floor((Y - miny)*scaley);
if (iy <0 || iy >= ny) return;
data[iy + ny*ix] += Weight;
}

double Histo2d::MaxBin(double &X, double &Y)
{
float *p, *pend;
int imax=0;
float valmax = -1e30;

for (p = data, pend = p + nx*ny; pend-p ; p++ )
  {
    if (*p > valmax) {valmax = *p; imax = p-data;}
  }
int ix = imax/ny;
int iy = imax - ix * ny;
X = minx + ((float)ix + 0.5)/scalex;
Y = miny + ((float)iy + 0.5)/scaley;
return valmax;
}

void Histo2d::ZeroBins(double &X, double &Y)
{
  int ix, iy;
  ix = (int) floor(( X - minx)*scalex);
  if (ix <0 || ix >= nx) return;
  iy = (int) floor((Y - miny)*scaley);
  if (iy <0 || iy >= ny) return;
  data[iy + ny*ix] = 0;
}
