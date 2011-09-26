#include "analyticpsf.h"

#include <vector>
#include <string>
using namespace std;

#include <cmath>


/*! 
  For regular users, there is no need to worry about what is going on 
  here: the ImagePSF class handles the needed interfaces.

  AnalyticPSF is a virtual class with actual derived classes that
  implement a given profile. You can choose your PSF type by name using
  ChooseAnalyticPSF. You can also add your own PSF to the list,
  without even altering this file. There is an example of test
  code (in order to make sure that your analytical derivatives
  are OK) at the end of this file. A trick is that the sign
  of the position derivatives has to be inverted, because the user
  worries about derivatives of profile(x_pix-x_star, ...) w.r.t
  x_star.

    The implemented derived classes have an integral of 1 over
  the whole plane. The constraint is integrated into the 
  derivatives w.r.t the profile parameters. This has to be enforced
  if you add new profiles.

  Analytic PSF's are integrated over pixels using Gauss quadrature,
 (as in DAOPHOT, where the relevant coefficients were taken),
  in the PixValue routine. The number of integration steps
  is choosen via a define (NPT)

*/


/* weights and abcissa for gauss integrations (borrowed from DAOPHOT),
   explained in Numerical recipes */

static double Dx[4][4] ={{0.00000000,  0.0,        0.0       , 0.0       },
			 {-0.28867513,  0.28867513, 0.0       , 0.0       },
			 {-0.38729833,  0.00000000, 0.38729833, 0.0       },
			 {-0.43056816, -0.16999052, 0.16999052, 0.43056816}};
static double Wt[4][4]= {{1.00000000,  0.0       , 0.0       , 0.0       },
			 {0.50000000,  0.50000000, 0.0       , 0.0       },
			 {0.27777778,  0.44444444, 0.27777778, 0.0       },
			 {0.17392742,  0.32607258, 0.32607258, 0.17392742}};
 




// number of points (per coordinate) to integrate over a pixel.
#define NPT 3


double AnalyticPSF::PixValue(const double &Xc, const double &Yc,
			     const double &XPix, const double &YPix,
			     const Vect &Params,
			     Vect *PosDer,
			     Vect *ParamDer) const
{
  double xPixCenter = floor(XPix+0.5);
  double yPixCenter = floor(YPix+0.5);
  Vect posDer(2);
  int npar = NPar();
  Vect paramDer(npar);
  double val = 0;
  for (int ix=0; ix<NPT; ++ix)
    {
      double x = xPixCenter+Dx[NPT-1][ix];
      double wx = Wt[NPT-1][ix];
      for (int iy=0; iy<NPT; ++iy)
	{
	  double y = yPixCenter+Dx[NPT-1][iy];
	  double weight = wx*Wt[NPT-1][iy];	  
	  val += weight*Profile(x-Xc,y-Yc, Params, PosDer, ParamDer);
	  if (PosDer)
	    {
	      posDer(0) += (*PosDer)(0)*weight;
	      posDer(1) += (*PosDer)(1)*weight;
	    }
	  if (ParamDer)
	    for (int ipar = 0; ipar<npar; ++ipar) 
	      paramDer(ipar) += weight*(*ParamDer)(ipar);
	}
    }
  if (PosDer) {(*PosDer)(0) = posDer(0); (*PosDer)(1) = posDer(1);}
  if (ParamDer) 
    for (int ipar = 0; ipar<npar; ++ipar) 
      (*ParamDer)(ipar) = paramDer(ipar);
#ifdef DEBUG
  cout << "PixValue " 
       << Xc << ' ' << Yc  << ' '
       << XPix << ' ' << YPix << ' ' 
       << Params
       << val << ' ' 
       << endl;
#endif
  return val;
}



static double sq(const double &x) { return x*x;};

/***************** GaussPSF *************************/

//! gaussian PSF with 3 parameters :wxx, wyy , wxy.
class GaussPSF : public AnalyticPSF {
 public :

  GaussPSF() {paramNames.push_back("wxx"); 
    paramNames.push_back("wyy");paramNames.push_back("wxy");}

  string Name() const {return "GAUSSIAN";}
  unsigned NPar() const { return 3;}
  double Profile(const double &X, const double &Y,
			 const Vect &Params,
			 Vect *PosDer = 0,
			 Vect *ParamGradient = 0) const;
  void InitParamsFromSeeing(const double &Seeing, Vect &Params) const;


  bool CheckParams(const Vect &Params) const
  { return (Params(0)*Params(1)-sq(Params(2))) > 0;}
		   


  virtual ~GaussPSF(){}; // warning killer

};
  


double GaussPSF::Profile(const double &X, const double &Y,
			  const Vect &Params,
			  Vect *PosDer,
			  Vect *ParamDer) const
{
  /* parameters (Params) wxx, wyy, wxy */
  double det = Params(0)*Params(1)-Params(2)*Params(2);
  double norm = sqrt(det)/(2*M_PI);
  double val = exp(-0.5*(X*X*Params(0)+Y*Y*Params(1)+2*X*Y*Params(2)))*norm;
  if (PosDer) // wrong sign on purpose (derivatives w.r.t -X)
    {
      (*PosDer)(0) = (X*Params(0)+Y*Params(2))*val;
      (*PosDer)(1) = (Y*Params(1)+X*Params(2))*val;
    }
  if (ParamDer)
    {
      (*ParamDer)(0) = -0.5*val*(X*X - Params(1)/det);
      (*ParamDer)(1) = -0.5*val*(Y*Y - Params(0)/det);
      (*ParamDer)(2) = -val*(X*Y  - Params(2)/det);
    }
  return val;
}
  
  
void GaussPSF::InitParamsFromSeeing(const double &Seeing, Vect &Params) const
{
  if (Params.Size() < NPar()) Params = Vect(NPar());
  Params(0) = 1/(sq(Seeing));
  Params(1) = Params(0);
  Params(2) = 0;
}


/************* MoffatPSF *******************/


class MoffatPSF : public AnalyticPSF {
protected :
  const double exponent;
 public :

  MoffatPSF(const double &Exponent) : exponent(Exponent) 
  {paramNames.push_back("wxx"); 
    paramNames.push_back("wyy");paramNames.push_back("wxy");}

  unsigned NPar() const { return 3;}
  double Profile(const double &X, const double &Y,
			 const Vect &Params,
			 Vect *PosDer = 0,
			 Vect *ParamGradient = 0) const;
  void InitParamsFromSeeing(const double &Seeing, Vect &Params) const;

  bool CheckParams(const Vect &Params) const
  { return (Params(0)*Params(1)-sq(Params(2))) > 0;}


  virtual ~MoffatPSF(){}; // warning killer

};


/**************** Moffat20PSF ****************/

class Moffat20PSF : public MoffatPSF {
 public :

  Moffat20PSF() :MoffatPSF(2.0) {};

  string Name() const {return "MOFFAT20";}

  virtual ~Moffat20PSF(){}; // warning killer

};

class Moffat25PSF : public MoffatPSF {
 public :

  Moffat25PSF() :MoffatPSF(2.5) {};

  string Name() const {return "MOFFAT25";}

  virtual ~Moffat25PSF(){}; // warning killer

};

class Moffat30PSF : public MoffatPSF {
 public :

  Moffat30PSF() :MoffatPSF(3.0) {};

  string Name() const {return "MOFFAT30";}

  virtual ~Moffat30PSF(){}; // warning killer

};
  


double MoffatPSF::Profile(const double &X, const double &Y,
			  const Vect &Params,
			  Vect *PosDer,
			  Vect *ParamDer) const
{
  /* parameters (Params) wxx, wyy, wxy */
  double det = Params(0)*Params(1)-Params(2)*Params(2);
  // in principle we should use sqrt(det). fabs saves some nan.
  double norm = sqrt(fabs(det))*(exponent-1)/M_PI;
  double fact =  1./(1.+(X*X*Params(0)+Y*Y*Params(1)+2*X*Y*Params(2)));
  double val = pow(fact,exponent)*norm;
  if (PosDer)
    {
      (*PosDer)(0) = 2*exponent*val*fact*(X*Params(0)+Y*Params(2));
      (*PosDer)(1) = 2*exponent*val*fact*(Y*Params(1)+X*Params(2));
    }
  if (ParamDer)
    {
      (*ParamDer)(0) = -exponent*val*fact*(X*X)+ 0.5*val*Params(1)/det;
      (*ParamDer)(1) = -exponent*val*fact*(Y*Y)+ 0.5*val*Params(0)/det;
      (*ParamDer)(2) = -2*exponent*val*fact*(X*Y) - val*Params(2)/det;
    }
  return val;
}

void MoffatPSF::InitParamsFromSeeing(const double &Seeing, Vect &Params) const
{
  if (Params.Size() < NPar()) Params = Vect(NPar());
  Params(0) = sq(0.6547/Seeing);
  Params(1) = Params(0);
  Params(2) = 0;
}


/***********************************************************************/


vector<const AnalyticPSF*> Profiles;


void AddAnalyticPSF(const AnalyticPSF *PSF)
{
  Profiles.push_back(PSF);
}




const AnalyticPSF *ChooseAnalyticPSF(const string &Name)
{
  if (Profiles.size() == 0)
    {
      AddAnalyticPSF(new GaussPSF());
      AddAnalyticPSF(new Moffat20PSF());
      AddAnalyticPSF(new Moffat25PSF());
      AddAnalyticPSF(new Moffat30PSF());
    }
  for (unsigned k=0; k < Profiles.size(); ++k)
    if (Profiles[k]->Name() == Name) return Profiles[k];
  cout << "ChooseAnalyticPSF : no profile called " <<  Name << endl;
  cout << " list of available profiles  :" << endl;
  for (unsigned k=0; k < Profiles.size(); ++k) 
    cout << Profiles[k]->Name() << endl; 
  cout << " ************** " << endl;
  return NULL;
}


#ifdef TEST_CODE

/* code to test a psf profile (in order to make sure that the 
derivatives are OK....)
*/


#include "analyticpsf.h"

int main( int nargs, char **args)
{
  if (nargs != 2)
    {
      cout << " usage : " << args[0] << " <PSFName> " << endl;
      exit(-1);
    }
  const AnalyticPSF *psf = ChooseAnalyticPSF(args[1]);
  unsigned npar = psf->NPar();

  cout << "# x : "<< endl << "# y: " << endl << "# v: " << endl 
       << "#dvdxexp : " << endl << "# dvdxth : " << endl
       << "#dvdyexp : " << endl << "# dvdyth : " << endl;
  for (unsigned k=0; k < npar; ++k)
    cout << "#grexp_" << k << " : " << endl << "#grth_" << k << " : " << endl;
  cout << "#end" << endl;




  double xc = 0.85,yc = 0.70;
  
  Vect posDer(2);
  Vect paramDer(npar);
  Vect params;
  psf->InitParamsFromSeeing(2, params);
  params(1) = params(0)*2; // introduce some asymmetry 
  params(2) = 0.1*params(0); // introduce some asymmetry 

  for (int i=-10; i <=10 ; ++i)
    for (int j=-10 ; j <= 10; ++j)
      {
	double v = psf->PixValue(xc,yc, i, j, params, &posDer, &paramDer);
	// posder
	double eps = 0.001;
	double vx = psf->PixValue(xc+eps, yc, i, j, params);
	double dvdxexp = (vx-v)/eps;
	double vy = psf->PixValue(xc, yc+eps, i, j, params);
	double dvdyexp = (vy-v)/eps;
	Vect grad(npar);
	for (unsigned k=0; k< npar; ++k)
	  {
	    Vect poff(params);
	    poff(k) += eps;
	    double vo = psf->PixValue(xc,yc, i, j, poff);
	    grad(k) = (vo-v)/eps;
	  }
	    
	cout << i << ' ' << j << ' '
	     << v << ' '
	     << dvdxexp << ' ' << posDer(0) << ' '
	     << dvdyexp << ' ' << posDer(1) << ' ';
	for (unsigned k=0; k < npar; ++k) 
	  cout << grad(k) << ' ' << paramDer(k) << ' ';
	  cout << endl;
      }

  cerr << " When analyzing the data, do not forget to check that the integral is reasonnably close to 1. " << endl;
 

  return EXIT_SUCCESS;
}

#endif /* TEST_CODE */

