#ifndef ANALYTICPSF__H
#define ANALYTICPSF__H


#include "matvect.h"

#include <string>
#include <vector>

using namespace std;





//! The virtual (interface) class. It holds the PixValue routine that integrates over pixels.
class AnalyticPSF {
 protected :
  //  ProfileFunc *profile;
  // const int npar;
  vector<string> paramNames;


 public :

  //     AnalyticPSF(const string &Name, ProfileFunc *Profile, 
  //      const int NPar , const vector<string> ParamNames);

  virtual string Name() const = 0;
  virtual string ParamName(const int Rank) const {return paramNames.at(unsigned(Rank));};
  virtual unsigned NPar() const  = 0;

 //! integrates PSF and requested derivatives over the pixel that contains XPix and YPix (pixel limits are at integer values + 1/2)
  double PixValue(const double &Xc, const double &Yc,
		  const double &XPix, const double &YPix,
		  const Vect &Params,
		  Vect *PosDer = 0,
		  Vect *ParamDer = 0) const;

  virtual double Profile(const double &X, const double &Y,
			 const Vect &Params,
			 Vect *PosDer = 0,
			 Vect *ParamGradient = 0) const = 0;

  virtual void InitParamsFromSeeing(const double &Seeing, Vect &Params) const = 0;

  virtual bool CheckParams(const Vect &Params) const = 0;

  virtual ~AnalyticPSF() {};
};


void AddAnalyticPSF(const AnalyticPSF *PSF);
const AnalyticPSF *ChooseAnalyticPSF(const string &Name);  

#endif /* ANALYTICPSF__H */
