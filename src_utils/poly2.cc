// -*- C++ -*-
#include <iostream>

#include "poly2.h"

using namespace std;

static double inline scal_prod(const double *x, const double *y, const unsigned n)
{
  register double val = 0;
  for (unsigned k=n; k; k--) { val += (*x) * (*y); ++x; ++y;}
  return val;
}



Poly2::Poly2(const double XMin, const double YMin, 
	const double XMax, const double YMax,
	     const int Deg) 
  : deg(Deg), nterms((deg+1)*(deg+2)/2), coeffs(nterms)
{

  // center and scale x and y ([-1,1]) to avoid numerical problems when fitting
  ax = (XMax-XMin)? 2./(XMax-XMin) : 1;
  bx = 1-ax*XMax;
  ay  = (YMax-YMin)? 2./(YMax-YMin) : 1;
  by = 1-ay*YMax;
}

void Poly2::DecreaseDegree()
{
  if (deg == -1) return;
  deg--;
  nterms =(deg+1)*(deg+2)/2;
} 

void Poly2::Monomials(const double &X, const double &Y, Vect &M) const
{
  int k=0;
  double x = ax*X+bx;
  double y = ay*Y+by;
  double xx = 1;
  for (int ix = 0; ix<=deg; ++ix)
    {
      double yy = 1;
      for (int iy = 0; iy<=deg-ix; ++iy)
	{
	M(k++) = xx*yy;
	yy *= y;
      }
    xx *= x;
    }
}

double Poly2::Value(const double &X, const double &Y) const
{
  Vect monom(nterms);
  Monomials(X, Y, monom);
  return scal_prod(monom.Data(), coeffs.Data(), nterms);
}


void Poly2::SetCoeffs(const Vect &Coeffs)
{
  if (Coeffs.Size() < nterms)
    {
      cerr << " Poly2::SetCoeffs : not enough coefficients provided (fatal)" 
	   <<  endl;
      abort();
    }
  // hand copy, in case Coeffs.Size > nterms.
  for (unsigned k=0; k < nterms; ++k) coeffs(k) = Coeffs(k);
}

void Poly2::SetCoeffs(const double *Coeffs)
{
  for (unsigned k=0; k < nterms; ++k) coeffs(k) = Coeffs[k];
}


void Poly2::Write(std::ostream &S) const
{
  S << "Poly2_version 1" << endl;
  S << ax << ' ' << bx << ' ' << ay << ' ' << by << endl;
  S << deg << endl;
  coeffs.writeASCII(S);
}


bool Poly2::Read(std::istream &S)
{
  string format;
  int version;
  S >> format >> version;
  if (format != "Poly2_version" || version != 1)
    {
      cout << " Poly2::Read : missing format line " << endl;
      return false;
    }
  S >> ax >> bx >> ay >> by;
  S >> deg;
  nterms = (deg+1)*(deg+2)/2;
  coeffs.readASCII(S);
  if (coeffs.Size() != nterms)
    {
      cout << " Poly2::Read : inconsistency between redundant data " << endl;
      return false;
    }
  return true;
}



std::ostream& operator << (std::ostream &S, const Poly2 &P)
{   P.Write(S); return S; }

std::istream& operator >> (std::istream &S, Poly2 &P)
{   P.Read(S); return S; }
