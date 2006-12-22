#ifndef POLY1__H
#define POLY1__H

#include <string>

#include "matvect.h"

#include "fstream"



//! Regular 1D polynomial.
class Poly1
{
 private :
  Vect params;
  unsigned npar;
  double a,b;

  /*
    a and b are chosen to map the x varying interval on [-1,1] (fit numerics)
    val(x) = sum_k params(k)*(a*x+b)**k
  */

 public:
  
 Poly1(const int Deg, const double XMin, const double XMax):
  npar(Deg+1)   
  {
    a = (XMax-XMin)? 2./(XMax-XMin) : 1;
    b = 1-a*XMax;
  };

  Poly1(const std::string &FileName) { npar = 0; a=b=0; Read(FileName);}

  void ParamDerivatives(const double &X, Vect &H) const
  {
    double xred = a*X+b;
    double p = 1;
    for (unsigned k=0; k <npar; ++k)
      {
	H(k) = p;
	p*= xred;
      }
  }

  double Value(const double &X) const
  {
    double val = 0;
    Vect der(npar);
    ParamDerivatives(X,der);
    for (unsigned k=0; k <npar; ++k) val += der(k)*params(k);
    return val;
  }

  unsigned NPar() const { return npar;};

  void SetParams(const Vect &P)
  {
    params.allocate(npar);
    for (unsigned k =0; k < npar; ++k) params(k) = P(k);
  }
  
  void Write(std::ostream &file) const;
  void Write(const std::string &FileName) const;

  const Vect& Params() const { return params;}

  void operator += (const Poly1 &Right);

  bool Read(const std::string &FileName);

  bool Read(std::ifstream &file, const std::string &FileName);




};




#endif /* POLY1__H */
