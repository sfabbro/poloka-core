#ifndef NONLINMODEL__H
#define NONLINMODEL__H

#include "poly1.h"
#include "string"

//! A class to describe non linearity of CCD response.
/*! the idea is simple: non linearity should preserve the value "0"
  and the derivative in 0. A regular polynomial time X^2 does exactly
  that.  Note that this class models the DIFFERENCE to non linearity,
  i.e. the Value() is 0 for a perfectly linear response.
*/
class NonLinModel : public Poly1
{
 private : 
  double maxPix;
  Mat cov;
  

 public :
  NonLinModel(const int Deg, const double XMin, const double XMax) : Poly1(Deg,XMin,XMax), maxPix(0), cov(NPar(), NPar()) {};


  NonLinModel(const std::string& FileName);

  double Value(const double &X) const { return (X*X*Poly1::Value(X));}

  void ParamDerivatives(const double &X, Vect &Der) const
  {
    Poly1::ParamDerivatives(X, Der);
    Der *= (X*X);
  }
  

  double MaxPix() const { return maxPix;}

  void SetMaxPix(const double &Val) { maxPix = Val;}


  void CumulateValues(const Poly1 &Right) { (Poly1 &)(*this) += Right; }

  bool Write(const std::string& FileName) const;

  bool Read(const std::string& FileName);

  const Mat& Cov() const { return cov;}

  double Cov(unsigned i, unsigned j) const { return cov(i,j);}
  
  Mat& Cov() {return  cov;}


};




#endif /* NONLINMODEL__H */
