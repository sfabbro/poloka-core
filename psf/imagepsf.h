#ifndef IMAGEPSF__H
#define IMAGEPSF__H


#include "reducedimage.h"
#include "countedref.h"

class PSFStarList;
class AnalyticPSF;
class Kernel;
class NonLinModel;

#include "matvect.h"

class Poly2;
class Poly2Image;


class ImagePSF : public RefCount {
 private :
  const ReducedImageRef reducedImage; /* imposes that reducedImage was allocated by new. Corresponding comment placed on constructor "documentation". */
  bool refitPSF;
  const AnalyticPSF* analyticPSF;
  int hSizeX, hSizeY;
  double seeing, gain, skylev;
  int nx,ny;
  int npar;
  vector<Poly2> psfParams;
  Poly2Image* residuals;
  NonLinModel* nonLinearity;
  

 public :
  
  //! Constructor. RI should have been allocated via "new".
  ImagePSF(const ReducedImage &RI, bool RefitPSF= false, const string & psffile_name="");

  //! half size of the PSF in pixels
  int HSizeX() const { return hSizeX;}

  //! half size of the PSF in pixels
  int HSizeY() const { return hSizeY;}

  //! limits of the stamp around Where StartI is in the stamp, EndI is beyond the last column.
  void StampLimits(const double &X, const double &Y,
		   int &BeginI, int &EndI,
		   int &BeginJ, int &EndJ) const;

  //! Seeing of the reduced image
  double Seeing() const { return seeing;}

  //! Gain of the reduced image
  double Gain() const { return gain;}

  //! number of parameters of the chosen analytic PSF
  unsigned NPar() const { return npar;}
  
  //! caries out the PSF fit
  bool FitPSF(PSFStarList &Stars);

  //! Access to the current PSF, with user provided Params.
  double PSFValueWithParams(const double &Xc, const double &Yc, 
			    const int IPix, const int JPix,
			    const Vect &Params,
			    Vect *PosDer, Vect *ParamDer, double *AnalyticValue = 0) const;

  //! Access to the current PSF pixels.
  double PSFValue(const double &Xc, const double &Yc, 
		  const int IPix, const int JPix,
		  Vect *PosDer = 0, Vect *ParamDer = 0, double *AnalyticValue = NULL) const;


  //! Access to current analytical PSF params (which may depend on position in the frame).
  Vect PSFParams(const double &X, const double &Y) const;

  //! residual image of fitted stars.
  void ResidualImage(const string &ResFitsName, const PSFStarList &List) const;

  // Residual tuple of fitted stars
  void ResidualTuple(const string &ResTupleName, 
		     const PSFStarList &Stars)const;

  ~ImagePSF();

  void Write(const std::string &FileName) const;

  void Write() const;

  bool FitNonLinearity(const bool WriteResTuple) /* const */;


 private :
  bool StackResiduals(PSFStarList &Stars, 
		      const Image &I, const Image &W, NonLinModel *NonLin = NULL);

  void ApplyNonLinearities(Image &I, const NonLinModel& NonLin) const;

  bool FitParametersVariation(PSFStarList &Stars, double MaxChi2);

  void Write(std::ostream &S) const;
  bool Read(std::istream &S);

  void Chi2Stat(const PSFStarList &Stars, double &mean, 
		double &median, double &sigma, 
		double &chi2_tot, int &ndof, const string &Message) const;

  int RemoveOutlierPixels(const PSFStarList &Stars,const Image &I, Image &W, const double NSigResCut) const;

  void SetStarsChi2(const Image &I, const Image &W, PSFStarList &Stars) const;

 public :
  void test_derivatives(const Image &I) const;

};

void SetPSFCardsFileName(const string &FileName);

bool MakePSF(const string &ImageName, const bool RefitPSF = false,
	     const bool UseExternalCatalog = true, const string & ExternalCatalogName="");

#endif /* IMAGEPSF__H */
