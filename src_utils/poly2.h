// -*- C++ -*-
#ifndef POLY2_H
#define POLY2_H

#include <fstream>
#include "matvect.h"



class Poly2 
{
 private:
  double ax, bx, ay, by;
  int deg;
  unsigned nterms;
  Vect coeffs;

 public:
  Poly2(const double XMin, const double YMin, 
	const double XMax, const double YMax,
	const int Deg);

  //! avoid using this one. assumes an immediatly following Read.
  Poly2() {};

  void Monomials(const double &X, const double &Y, Vect &M) const;
  double Value(const double &X, const double &Y) const;
  
  void SetCoeffs(const Vect &Coeffs);
  void SetCoeffs(const double *Coeffs);
  unsigned NTerms() const { return nterms;}

  bool Read(std::istream &S);
  void Write(std::ostream &S) const;

  void DecreaseDegree();

};


std::ostream& operator << (std::ostream &S, const Poly2 &P);
std::istream& operator >> (std::istream &S, Poly2 &P);



#endif
