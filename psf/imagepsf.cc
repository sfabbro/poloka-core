#include "imagepsf.h"
#include "analyticpsf.h"
#include "psfstar.h"
#include "vutils.h"
#include "wcsutils.h" // for WCSFromHeader
#include "gtransfo.h"
#include "fileutils.h"

#include "dimage.h" // for Kernel
#include "poly2.h"
#include "poly1.h" // used for non linearity of response fitting.

#include "polokaexception.h"

#include <cmath> //for floor(double)

/*! Overall structure of the PSF modelling:

The heart of the system is the class ImagePSF, where the driver
routine is FitPSF. There are some possibilities to provide parameters
and get more output via datacards.

  The PSF model is 2 folds: an analytic PSF, and some discretized
residuals. Both can be variable in the field (using polynomials
of choosable degrees). There are basically 2 modes, one to
fit the PSF, and one to use it for fits.


- auxillary classes:

 PSFCards: data cards for this module.
 Poly2 : handler for 2 dimensional polynomials.
 Poly2Image : same as above, but for an image where each pixel
 is itself a polynomial.
 NonLinModel : to describe non-linearity of CCD response.


 AnalyticPSF : see analyticpsf.cc for comments

 PSFStar(List) : really auxillary. it only holds PSF specific parameters
 no io's beyond what is needed for a control tuple.

*/


/* TODO
   - degres de variation spatiales bornees par la longueur de la liste
   - swap indices in filling loop of StackResiduals 

*/

#define PSF_FILE_NAME "psf.dat"


struct PSFCards
{
  /* when you add a field here, your should take care of it in 3 other places:
     - ReadCards
     - dump
     - PSFCards::PSFCards
     
     A more centralized handling would be welcome (as e.g. in SExtractor)
  */

  double nSigPSFChi2Cut;
  double nSigPSFHalfSize;
  double nSigPixResCut;
  //  int maxPSFHalfSize;
  int spatialAnalyticVarDeg;
  int spatialPixVarDeg;
  int nonLinearityDegree;
  string analyticKind;
  string outputDirectory; 
  string allResidualsFileName;
  string lastResidualsFileName;
  string residualsTupleFileName;
  string residualsImageFileName;

#define WRITE_IF_DEFINED(STREAM,NAME,VAR) \
  if(VAR!="") STREAM << NAME << VAR << endl;

  PSFCards();
  bool ReadCards(const string&);
  void dump(ostream &s = cout)
  {
    s << " PsfCards : " << endl;
    s << "PSF_SPATIAL_ANALYTIC_VAR_DEG " << spatialAnalyticVarDeg << endl;
    s << "PSF_SPATIAL_PIX_VAR_DEG " << spatialPixVarDeg << endl;
    s << "PSF_NSIG_CHI2_CUT " << nSigPSFChi2Cut << endl;
    s << "PSF_NSIG_HALF_SIZE " << nSigPSFHalfSize << endl;
    s << "PSF_NSIG_PIX_RES_CUT " << nSigPixResCut << endl;
    //    s << "MAX_PSF_HALF_SIZE" << maxPSFHalfSize << endl;
    s << "PSF_ANALYTIC_KIND " << analyticKind << endl; 
    s << "PSF_NON_LINEARITY_DEGREE " << nonLinearityDegree << endl; 
    WRITE_IF_DEFINED(s,"PSF_OUTPUT_DIR ",outputDirectory);
    WRITE_IF_DEFINED(s, "PSF_ALL_RESIDUALS_FILE_NAME ",allResidualsFileName);
    WRITE_IF_DEFINED(s,"PSF_LAST_RESIDUAL_FILENAME ",lastResidualsFileName);
    WRITE_IF_DEFINED(s,"PSF_RESIDUALS_TUPLE_FILENAME ",residualsTupleFileName);
    WRITE_IF_DEFINED(s,"PSF_RESIDUALS_IMAGE_FILENAME ",residualsImageFileName);
  }
};

#define READ_IF_EXISTS(VAR,TAG,TYPE) \
if (cards.HasKey(TAG)) VAR=cards.TYPE(TAG)

#include "datacards.h"
#include "toadscards.h"

static string PSFCardsFileName;

void SetPSFCardsFileName(const string &FileName)
{
  PSFCardsFileName = FileName;
}


bool PSFCards::ReadCards(const string &FileName)
{
  if (!FileExists(FileName))
    {
      std::cerr << " cannot open datacards" << FileName << std::endl;
      return false;
    }  
  DataCards cards(FileName);
  cout << "Reading PSF datacards in " << FileName << endl;
  READ_IF_EXISTS(spatialAnalyticVarDeg,"PSF_SPATIAL_ANALYTIC_VAR_DEG",IParam);
  READ_IF_EXISTS(spatialPixVarDeg,"PSF_SPATIAL_PIX_VAR_DEG",IParam);
  READ_IF_EXISTS(nSigPSFChi2Cut,"PSF_NSIG_CHI2_CUT",DParam);
  READ_IF_EXISTS(nSigPSFHalfSize,"PSF_NSIG_HALF_SIZE",DParam);
  READ_IF_EXISTS(nSigPixResCut,"PSF_NSIG_PIX_RES_CUT",DParam);
  READ_IF_EXISTS(nonLinearityDegree,"PSF_NON_LINEARITY_DEGREE",IParam);

  //  READ_IF_EXISTS(maxPSFHalfSize,"MAX_PSF_HALF_SIZE",IParam);
  READ_IF_EXISTS(analyticKind,"PSF_ANALYTIC_KIND",SParam);
  READ_IF_EXISTS(outputDirectory,"PSF_OUTPUT_DIR",SParam);
  READ_IF_EXISTS(outputDirectory, "PSF_OUTPUT_DIR",SParam);
  READ_IF_EXISTS(allResidualsFileName,"PSF_ALL_RESIDUALS_FILE_NAME",SParam);
  READ_IF_EXISTS(lastResidualsFileName,"PSF_LAST_RESIDUAL_FILENAME",SParam);
  READ_IF_EXISTS(residualsTupleFileName,"PSF_RESIDUALS_TUPLE_FILENAME",SParam);
  READ_IF_EXISTS(residualsImageFileName,"PSF_RESIDUALS_IMAGE_FILENAME",SParam);
  return true;
}



PSFCards::PSFCards()
{
  spatialAnalyticVarDeg = 2;
  spatialPixVarDeg = 1;
  nSigPSFChi2Cut = 4;
  nSigPSFHalfSize = 5;
  nSigPixResCut = 5;
  nonLinearityDegree = 1;
  //  maxPSFHalfSize = 15;
  analyticKind = "MOFFAT25";
  if (PSFCardsFileName != "") ReadCards(PSFCardsFileName);
  else ReadCards(DefaultDatacards());
  dump();
}

/* Singleton management */
static PSFCards *ThePSFCards = NULL;

const PSFCards &Cards()
{
  if (!ThePSFCards) ThePSFCards = new PSFCards();
  return *ThePSFCards;    
}


#include "nonlinmodel.h"



/******************* handy local stuff ***************/

static double inline scal_prod(const double *x, const double *y, const unsigned n)
{
  register double val = 0;
  for (unsigned k=n; k; k--) { val += (*x) * (*y); ++x; ++y;}
  return val;
}

static double sq(const double &x) { return x*x;}



/************** ResampCoeffs ********************/


// Image Resampling Coefficients (no oversampling)
/* This provides coefficients along one coordinate. Just get 2 sets of
coefficients and multiply them to get actual image resampling.  We may
consider providing several image resampling kernels, but what for?
The key reason to extract it form the code flow is that we need it a 2
different places, and coding it only once ensures that it is the same
!
*/


// Bilinear is enough for what we need here
static int ResampNCoeffs() { return 2;}

static void ResampCoeffs(const double Dx, int &StartI, int &EndI, Vect &Coeffs,
			 Vect *DCoeffs = 0)
{
  if (Dx>=0)
    {
      StartI = 0;
      EndI = 1;
      Coeffs(0) = (1-Dx);
      Coeffs(1) = Dx;
    }
  else
    {
      StartI = -1;
      EndI = 0;
      Coeffs(0) = -Dx;
      Coeffs(1) = (1+Dx);
    }
  if (DCoeffs) { (*DCoeffs)(0) = -1; (*DCoeffs)(1) = 1;}
#ifdef STORAGE
  //nearest neighbor
  StartI = 0;
  EndI=0;
  Coeffs(0) = 1;
#endif
}    


/************** class Poly2Image **************************/
/* 
   a class to handle image pixels (here residuals of the PSF)
   which depend on 2 parameters (here coordinates in the frame)
*/

//#include "dimage.h"
//#include <vector>

class Poly2Image 
{
 private:
  double ax, bx, ay, by;
  int deg;
  unsigned nterms;
  int hSizeX, hSizeY;
  vector<Kernel> coeffs;
  Vect pixSum; // integrals of the kernels that hold the coefficients
  Vect biasX; // sum_ij(kernel(i,j)*i)
  Vect biasY; // sum_ij(kernel(i,j)*j)


 public:
  
  Poly2Image(const double XMin, const double YMin, 
	     const double XMax, const double YMax,
	     const int Deg, const int HSizeX, const int HSizeY);

  // don't use this : it assumes an immediately following Read.
  Poly2Image( const int HSizeX, const int HSizeY) 
    : hSizeX(HSizeX), hSizeY(HSizeY) {};

  void Monomials(const double &X, const double &Y, Vect &M) const;
  double PixValue(const double &X, const double &Y, 
		  const int I, const int J, 
		  Vect *PosDer = 0) const;

  //! integral and "unormalized averages" (likely to be useless) at position Xc Yc
  void Moments(const double &Xc, const double &Yc, double &Integral, double &Bx, double &By) const;

  // Kernel Value(const double &X, const double &Y) const;

  //! assumes that the ordering of pixels is the one of "Kernel".
  void SetCoeffs(const Vect &Coeffs);

  unsigned NTerms() const { return nterms;}
  int NPix() const { return (2*hSizeX+1)*(2*hSizeY+1);}

  void operator += (const Poly2Image &R);
  const Kernel& Coeffs(const int Rank) const { return coeffs[Rank];}

  void Write(std::ostream &S) const;
  bool Read(std::istream &S);

#ifdef STORAGE
  typedef enum ConstraintType {Integral, Position, IntegralAndPosition};

  void Poly2Image::ComputeConstraintMat(Poly2Image::ConstraintType Type,
					Mat &Constraint);

#endif

};

Poly2Image::Poly2Image(const double XMin, const double YMin, 
	const double XMax, const double YMax,
	     const int Deg, const int HSizeX, const int HSizeY)
  : deg(Deg), nterms((deg+1)*(deg+2)/2), hSizeX(HSizeX), hSizeY(HSizeY)
{

  // center and scale x and y ([-1,1]) to avoid numerical problems when fitting
  ax = (XMax-XMin)? 2./(XMax-XMin) : 1;
  bx = 1-ax*XMax;
  ay  = (YMax-YMin)? 2./(YMax-YMin) : 1;
  by = 1-ay*YMax;
}

void Poly2Image::Monomials(const double &X, const double &Y, Vect &M) const
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


void Poly2Image::Moments(const double &Xc, const double &Yc, 
			 double &Integral, double &Bx, double &By) const
{
  Vect monom(nterms);
  Monomials(Xc,Yc,monom);
  Integral = scal_prod(monom.Data(), pixSum.Data(), nterms);
  Bx = scal_prod(monom.Data(), biasX.Data(), nterms);
  By = scal_prod(monom.Data(), biasY.Data(), nterms);
}
  

#ifdef STORAGE
/* This routine is correct but useless: one should not apply
   this kind of constraints for a precise modelling of the PSF.
*/
void Poly2Image::ComputeConstraintMat(Poly2Image::ConstraintType WhichConstraint,
				      Mat &Constraint)
{
  unsigned nc=0;
  switch (WhichConstraint)
    {
    case Integral : nc = nterms; break;
    case Position : nc = 2*nterms; break;
    case IntegralAndPosition : nc = 3*nterms; break;
    }
      
  unsigned npix = NPix();
  Constraint.allocate(npix*nterms, nc);
  int ic = 0;

  if (WhichConstraint == Integral || WhichConstraint == IntegralAndPosition)
    {
      /* sum over pixels does not depend on position:
	 sum_i,j coeffs[k](i,j) = 0 for all k (but perhaps for k =0)
      */	 
      for (unsigned int k=0; k < nterms ; ++k)
	for (unsigned ipix = 0; ipix < npix; ++ipix)
	  Constraint(ipix + npix*k, k) = 1.;
      ic = nterms;
    }
  if (WhichConstraint == Position || WhichConstraint == IntegralAndPosition)
    {
      /* centroid does not depend on position:
	 sum_i,j i*coeffs[k](i,j) = 0, for all k 
	 sum_i,j j*coeffs[k](i,j) = 0, for all k 
      */
      for (unsigned int k=0; k < nterms ; ++k)
	{
	  unsigned c = 0;
	  for (int j=-hSizeY; j <= hSizeY; ++j)
	    for (int i=-hSizeX; i <= hSizeX; ++i)
	      {
		Constraint(c + k*npix, ic+k) = i; 
		Constraint(c + k*npix, ic+k+nterms) = j;
		++c;
	      }
	}
    }
} 
#endif /* STORAGE */


static void UnormalizedKernelMoments(const Kernel &K, double &Mx, double &My)
{
  Mx = 0;
  My = 0;
  for (int j=-K.HSizeY(); j <= K.HSizeY(); ++j)
  for (int i=-K.HSizeX(); i <= K.HSizeX(); ++i)
    {
      Mx += K(i,j)*double(i);
      My += K(i,j)*double(j);
    }
}


void Poly2Image::SetCoeffs(const Vect &Coeffs)
{
  int npix = NPix();
  if (Coeffs.Size() < npix*nterms)
    {
      cout << " Poly2Image::SetCoeffs : vector provided is too short " << endl;
      abort();
    }
  coeffs.clear();
  pixSum.allocate(nterms);
  biasX.allocate(nterms);
  biasY.allocate(nterms);
  for (unsigned k=0; k < nterms; ++k)
    {
      coeffs.push_back(Kernel(hSizeX, hSizeY));
      memcpy(coeffs.back().begin(), Coeffs.Data()+(k*npix), 
	     sizeof(double)*npix);
      pixSum(k) = coeffs.back().sum();
      UnormalizedKernelMoments(coeffs.back(), biasX(k), biasY(k));
    }
}

void Poly2Image::operator += (const Poly2Image &R)
{
  // non exhaustive compatibility check. if sizes disagree, Kernel::operator +=() will hopefully crash
  if ( fabs(ax-R.ax)>1e-5*ax  || nterms != R.nterms)
    {
      cout << "  Poly2Image::operator += : inconsistent definitions " << endl;
      abort();
    }
  pixSum += R.pixSum;
  for (unsigned k = 0; k < nterms; ++k)
    {
      coeffs[k] += R.coeffs[k];
      UnormalizedKernelMoments(coeffs[k], biasX(k), biasY(k));
      //DEBUG
      //      cout << " in Poly2Image::operator +=() sum bx by " << pixSum(k) << ' ' << biasX(k) << ' ' << biasY(k) << endl;
    }
}


/* There is some inefficiency associated to this routine since it is
 called once per pixel and most of the time (if not always), we want
 the whole stamp, and there are many calculations that hold for the
 whole stamp. However, this routine is 3 times faster (per stamp) than
 most of the analytical PSF's. So going to stamps instead of pixel
 calculations would not speed-up much, I guess.
*/

//static int debug_flag = 0;

double Poly2Image::PixValue(const double &Xc, const double &Yc, 
			    const int IPix, const int JPix, 
			    Vect *PosDer) const
{
  Vect monom(nterms);
  Monomials(Xc, Yc, monom);
  Vect cx(ResampNCoeffs());
  Vect cy(ResampNCoeffs());
  Vect dcx(ResampNCoeffs());
  Vect dcy(ResampNCoeffs());

  int ic = int(floor(Xc+0.5));
  int jc = int(floor(Yc+0.5));
  int minDi, maxDi, minDj, maxDj;
  ResampCoeffs(ic-Xc,minDi,maxDi,cx, &dcx);
  ResampCoeffs(jc-Yc,minDj,maxDj,cy, &dcy);
  // up to here IPix and JPix have not been used.
  double val = 0;
  double valdx = 0;
  double valdy = 0;
  for (unsigned int k=0; k<nterms ; ++k) 
    {
      double val1 = 0;
      double val1dx = 0;
      double val1dy = 0;
      const Kernel &pixVals = coeffs[k];
      for (int dj = minDj; dj<= maxDj; ++dj)
	{
	  int j = JPix-jc+dj;      
	  /* handling of the sides : assume that the pixels beyond the
	     frame are equal to 0, as in StackResiduals
	  */
	  if (abs(j) > hSizeY) continue;
	  for (int di = minDi; di<= maxDi; ++di)
	    {
	      int i = IPix-ic+di;
	      if (abs(i) > hSizeX) continue;
	      val1 += cx(di-minDi)*cy(dj-minDj)*pixVals(i,j);
	      /* Code used to check that derivatives are correct
	      if (debug_flag >10)
		cout << " Poly2Image::PixValue IPix JPix i j cx cy " 
		     << IPix << ' '  << JPix << ' ' << i << ' ' << j << ' ' << cx(di-minDi) << ' ' << cy(dj-minDj) << endl;
	      */
	      val1dx += dcx(di-minDi)*cy(dj-minDj)*pixVals(i,j);
	      val1dy += cx(di-minDi)*dcy(dj-minDj)*pixVals(i,j);
	    }
	}
      val += val1*monom(k);
      valdx += val1dx*monom(k);
      valdy += val1dy*monom(k);
    }
  if (PosDer)
    {
      /* - sign?? yes because PosDer is the derivative w.r.t to Xc which
	 enters with a - sign in ResampCoeffs. Same trick in AnalyticPSF
      */
      (*PosDer)(0) -= valdx;
      (*PosDer)(1) -= valdy;
    }
  return val;
}

void static kernel_write(const Kernel &K, std::ostream &S)
{
  S << K.HSizeX() << ' ' << K.HSizeY() << endl;
  const double *pend = K.begin()+K.Nx()*K.Ny();
  for (const double *p = K.begin(); p <pend; ++p)
    S << *p << ' ';
  S << endl;
}

void static kernel_read(Kernel &K, std::istream& S)
{
  int hx,hy;
  S >> hx >> hy;
  K.Allocate(2*hx+1, 2*hy+1);
  const double *pend = K.begin()+K.Nx()*K.Ny();
  for (double *p = K.begin(); p < pend; ++p)
    S >> *p;
} 

static void  kernel_crop(Kernel &K, const int hSizeX, const int hSizeY)
{
  if (K.HSizeX() == hSizeX  && K.HSizeY() == hSizeY) return;
  Kernel copy(K);
  K.Allocate(2*hSizeX+1, 2*hSizeY+1);
  for (int j = -hSizeY; j <= hSizeY; ++j)
    for (int i = -hSizeX; i <= hSizeX; ++i)
    {
      K(i,j) = copy(i,j);
    }
}
  

void Poly2Image::Write(std::ostream &S) const
{
  S << "Poly2Image_version 1 " << endl;
  S << ax << ' ' << bx << ' ' << ay << ' ' << by << endl;
  S << deg << ' ' << hSizeX << ' ' << hSizeY << endl;
  S << coeffs.size() << endl;
  for (unsigned k=0 ; k < coeffs.size(); ++k)
    kernel_write(coeffs[k], S);
}

bool Poly2Image::Read(std::istream &S)
{
  string comment;
  int version = 0;
  S >> comment >> version;
  if (comment != "Poly2Image_version" || version != 1)
    {
      cout << " error reading Poly2Image (version line not found) " << endl;
      return false;
    }
  
  S >> ax >> bx >> ay >> by;
  double hx,hy;
  S >> deg >> hx >> hy;
  nterms = (deg+1)*(deg+2)/2;
  unsigned nkern;
  S >> nkern;
  coeffs.clear();
  pixSum.allocate(nkern);
  biasX.allocate(nkern);
  biasY.allocate(nkern);
  for (unsigned k=0 ; k < nkern; ++k)
    {
      coeffs.push_back(Kernel());
      kernel_read(coeffs.back(), S);
      /* 
	 This "croping" is due to the fact that some PSF's were
	 produced with Poly2Image's with essentially undefined
	 "borders": the values beyond the hSizeX, hSizeY of the PSF
	 itself were extremely noisy.  The idea is to correct that
	 here (plus Poly2Image::PixValue and StackResiduals
      */
      kernel_crop(coeffs.back(), hSizeX, hSizeY);
      pixSum(k) = coeffs.back().sum();
      UnormalizedKernelMoments(coeffs[k], biasX(k), biasY(k));
      //DEBUG
      //      cout << " in Poly2Image::Read sum bx by " << pixSum(k) << ' ' << biasX(k) << ' ' << biasY(k) << endl;
      
    }
  return true;
}


/*********************   ImagePSF *************************/

static string PSFFileName(const ReducedImage &RI)
{
  return AddSlash(RI.Dir())+PSF_FILE_NAME;
}


#include <fstream>

ImagePSF::ImagePSF(const ReducedImage &RI, bool RefitPSF) 
  : reducedImage(&RI), refitPSF(RefitPSF)
{
  analyticPSF = NULL;
  residuals = NULL;
  nonLinearity = NULL;
  skylev = 0;
  string psfFileName = PSFFileName(*reducedImage);
  if (FileExists(psfFileName) && !refitPSF)
    {
      ifstream s(psfFileName.c_str());
      cout << "Reading PSF in file : " << psfFileName.c_str() << endl;
      Read(s);
      s.close();
      return;
    }
  // if we get here, either there is no PSF around or refitPSF = true
  // if there is not PSF around, we'd better fit the psf, so
  refitPSF = true; 
  seeing = reducedImage->GFSeeing();
  if (seeing <= 0) 
    seeing = reducedImage->Seeing();
  gain = reducedImage->Gain();
  if (gain <= 0)
    throw (PolokaException(" for image "+reducedImage->Name()+", ImagePsf cannot work with negative or null gain"));
  nx = reducedImage->XSize();
  ny = reducedImage->YSize();
  hSizeX = int(Cards().nSigPSFHalfSize*seeing)+1;
  //  hSizeX = min(hSizeX, Cards().maxPSFHalfSize());
  hSizeY = hSizeX;
  npar = 0;

}

void ImagePSF::Write() const
{
  string psfFileName = AddSlash(reducedImage->Dir())+PSF_FILE_NAME;
  Write(psfFileName);
  if (nonLinearity)
    nonLinearity->Write((reducedImage->Dir()+"/nonlinpar.dat").c_str());
}

void ImagePSF::Write(const std::string &FileName) const
{
  ofstream s(FileName.c_str());
  cout << " Writing PSF file : " << FileName << endl;
  Write(s);
  s.close();
}

void ImagePSF::Write(std::ostream &S) const
{
  //  ReducedImage reducedImage; ??
  S << "ImagePSF_version 1" << endl;
  S << analyticPSF->Name() << endl;
  S << hSizeX << ' ' <<  hSizeY << endl;
  S << seeing << ' ' << gain << endl;
  S << nx << ' ' << ny << endl;
  S << psfParams.size() << endl;
  for (unsigned k = 0; k < psfParams.size(); ++k)
    S << psfParams[k];
  if (residuals) 
    {
      S << "pixellized_residuals" << endl;
      residuals->Write(S);
    }
  else
    S << "No pixellized residuals " << endl;
}
 


bool ImagePSF::Read(std::istream &S)
{
  string info;//  ReducedImage reducedImage; ??
  int version;
  S >> info >> version;
  if (info != "ImagePSF_version" || version != 1) 
    {
      cout << "ImagePSF::Read : bad format : " << info << ' ' << version << endl;
      return false;
    }
  S >> info;
  if (!(analyticPSF = ChooseAnalyticPSF(info)) || analyticPSF->Name() != info)
    {
      cout << " could not get the PSF profile back. deep shit " << endl;
      return false;
    }
  npar = analyticPSF->NPar();
  S >> hSizeX >> hSizeY;
  S >> seeing >> gain;
  S >> nx >> ny;
  unsigned count;
  S >> count;
  psfParams.clear();
  for (unsigned k = 0; k < count; ++k)
    {
      psfParams.push_back(Poly2());
      S >> psfParams.back();
    }
  if (residuals) 
    { delete residuals; residuals = NULL;}
  S >> info;
  if (info == "pixellized_residuals")
    {
      residuals = new Poly2Image(hSizeX, hSizeY);
      residuals->Read(S);
    }
  return true;
}


ImagePSF::~ImagePSF()
{
  if (residuals) delete residuals;
}


/*! for the moment, just forward the call to PSFValue to the
corresponding AnalyticPSF routine.  we'll add pixelized corrections
here if needed
*/
double ImagePSF::PSFValueWithParams(const double &Xc, const double &Yc, 
				    const int IPix, const int JPix,
				    const Vect &Params,
				    Vect *PosDer, Vect *ParamDer) const
{
  return analyticPSF->PixValue(Xc,Yc,IPix, JPix, Params, PosDer, ParamDer);
}


Vect ImagePSF::PSFParams(const double &Xc, const double &Yc) const
{
  Vect params(npar);
  for (int k =0; k < npar; ++k)
    params(k) = psfParams[k].Value(Xc,Yc);
  return params;
}

double ImagePSF::PSFValue(const double &Xc, const double &Yc, 
			  const int IPix, const int JPix, 
			  Vect *PosDer, Vect *ParamDer, double *AnalyticValue) const
{

  double val = 0;
  double residuals_integral = 0;
  /* the policy is to provide user with a centered PSF, i.e.
     with a PSF which has a centroid equal to the requested coordinate.
     dx and dy are the offsets necessary to acheive that. They are 
     NOT applied to address the "residual" part, because they are not
     used in the fit (in StackResiduals). If we put them in here,
     then they also have to be accounted for there.*/
  double dx = 0;
  double dy = 0;
  if (residuals)
    {
      residuals->Moments(Xc,Yc,residuals_integral, dx, dy);
      // no offsets (dx,dy) applied, see comment above.
      val += residuals->PixValue(Xc,Yc,IPix, JPix, PosDer);
    }


  Vect params(npar);
  for (int k =0; k < npar; ++k)
    params(k) = psfParams[k].Value(Xc,Yc);
  double anal_val = analyticPSF->PixValue(Xc-dx,Yc-dy,
					  IPix, JPix, params, 
					  PosDer, ParamDer);
  val += (1.-residuals_integral) * anal_val;
  if (AnalyticValue) *AnalyticValue = anal_val;

  return val;
}  

void ImagePSF::StampLimits(const double &X, const double &Y, 
			   int &BeginI, int &EndI,
			   int &BeginJ, int &EndJ) const
{
  int iPix = int(floor(X+0.5));
  int jPix = int(floor(Y+0.5));
  BeginI = max(iPix-hSizeX,0);
  BeginJ = max(jPix-hSizeY,0);
  EndI = min(iPix+hSizeX+1,nx);
  EndJ = min(jPix+hSizeY+1,ny);
}			   


static int erase_outliers(PSFStarList &Stars, double MaxChi2)
{
  // discard outliers
  int erasedCount = 0;
  for (PSFStarIterator i = Stars.begin(); i != Stars.end();)
    {
      PSFStar &s = **i;
      if (s.PSFChi2() > MaxChi2) 
	{
	  cout << " erasing outlier star at " << s.x << ' ' << s.y 
	       << " fmax " << s.Fluxmax() 
	       << " chi2 = " << s.PSFChi2() << endl;
	  i = Stars.erase(i); 
	  erasedCount++;
	}
      else ++i;
    }
  return erasedCount;
}



bool ImagePSF::FitParametersVariation(PSFStarList &Stars, double MaxChi2)
{
  int erasedcount ;
  double mean, median, sigma, chi2_tot;
  int ndof;
  cout << " now Fit Parameters Variation " << endl;
  do 
    {
      unsigned nstars = Stars.size();
      Poly2 spatialVar(0,0,nx,ny, 
		       Cards().spatialAnalyticVarDeg);
      // adjust spatial variation to number of stars
      while (spatialVar.NTerms() > nstars*2) 
	{
	  spatialVar.DecreaseDegree();
	  if (spatialVar.NTerms() == 0)
	    return false;
	}
      unsigned nterms = spatialVar.NTerms();
      Vect b(nterms*npar);
      Mat a(nterms*npar,nterms*npar);
      Vect h(nterms);
      for (PSFStarCIterator i = Stars.begin(); i != Stars.end(); ++i)
	{
	  const PSFStar &s = **i;
	  if (s.PSFChi2()> MaxChi2) continue;
	  spatialVar.Monomials(s.x, s.y, h);
	  const Mat &w = s.PSFParamsWeight();
	  Vect wp =  (s.PSFParamsWeight()*s.PSFParams());
	  
	  for (int ip1 = 0; ip1 < npar; ++ip1)
	    for (unsigned ih1 = 0; ih1< nterms; ++ih1)
	      {
		int ind1 = ip1*nterms+ih1;
		b(ind1) += h(ih1)*wp(ip1);	  
		for (int ip2 = ip1; ip2 < npar; ++ip2)
		  for (unsigned ih2 = 0; ih2< nterms; ++ih2)
		    a(ind1, ip2*nterms+ih2) += w(ip1,ip2)*h(ih1)*h(ih2);
	      }
	}
      
      if (cholesky_solve(a, b, "U")!= 0)
	{
	  cout << " ImagePSF::FitParametersVariation : could not solve for PSF spatial variation " 
	       << endl;
	  return false;
	}
      
      psfParams = vector<Poly2>(npar, spatialVar);
      for (int k=0; k < npar; ++k) 
	{
	  psfParams[k].SetCoeffs(&b(nterms*k));
	}
      
      // Set chi'2 to detect outliers
      for (PSFStarIterator i = Stars.begin(); i != Stars.end(); ++i)
	{
	  PSFStar &s = **i;
	  Vect delta = PSFParams(s.x, s.y) - s.PSFParams();
	  double chi2 = delta*(s.PSFParamsWeight()*delta);
	  s.SetPSFChi2(chi2);
	}
      Chi2Stat(Stars, mean, median, sigma, chi2_tot, ndof, "");
      // overwrite ndof because in this fit we did not use the pixels as 
      // Chi2Stat assumes.
      ndof = nstars*npar - b.Size();
      cout << " count , average chi2 (Spat. Var)" 
	   << Stars.size() << ' ' 
	   << mean << ' '
	   << endl;
      cout << " *** chi2/ndof " << chi2_tot << '/' << ndof << "=" 
	   << chi2_tot/ndof << endl;

      double maxChi2 = median + Cards().nSigPSFChi2Cut * sigma;
      erasedcount = erase_outliers(Stars,maxChi2);
    } while (erasedcount != 0);
  // output for grep in logs:
  cout << " PSF_SPATVAR_name_nstar_chi2_ndof_chi2/ndof "
       << reducedImage->Name() << ' '
       << Stars.size() << ' '
       << chi2_tot << ' '
       << ndof << ' '
       << chi2_tot/ndof << ' '
       << endl;
  return true;
}

void ImagePSF::Chi2Stat(const PSFStarList &Stars, double &mean, 
			double &median, double &sigma, 
			double &chi2_tot, int &ndof, 
			const string &Message) const
{
  unsigned nstars = Stars.size();
  double *chi2 = new double[nstars];
  int count = 0;
  for (PSFStarCIterator i = Stars.begin(); i != Stars.end(); ++i)
    chi2[count++] = (*i)->PSFChi2();
  Dmean_median_sigma(chi2, count, mean, median, sigma);
  // 3 is for x,y, flux
  ndof = nstars*((2*hSizeX+1)*(2*hSizeY+1) - 3 - analyticPSF->NPar());
  if (residuals)
    ndof -= residuals->NTerms() * residuals->NPix();
  chi2_tot = mean*nstars;
  delete [] chi2;
  if (Message != "")
    {
      cout << " count , average chi2 (" << Message << ") " 
	   << Stars.size() << ' ' 
	   << mean << ' '
	   << endl;
      cout << " *** chi2/ndof " << chi2_tot << '/' << ndof << "=" 
	   << chi2_tot/ndof << endl;
    }
}


double ImagePSF::ComputeChi2(const PSFStarList &Stars, const Image &I, const Image &W) const
{
  double chi2 = 0;
  for (PSFStarCIterator it = Stars.begin(); it != Stars.end(); ++it)
    {
      const PSFStar &s = **it;
      double xc = s.PSFX();
      double yc = s.PSFY();
      double psfFlux = s.PSFFlux();
      int starti, endi, startj,endj;
      StampLimits(xc,yc,starti, endi, startj, endj);
      for (int j=startj; j < endj; ++j)
	for (int i=starti; i < endi; ++i)
	  {
	    double w = W(i,j);
	    if (w<=0) continue;	    
	    double psfVal = PSFValue(xc,yc,i,j);
	    w = 1./(1./w+psfVal*psfFlux/gain);
	    double res = (I(i,j)-psfFlux*psfVal);
	    chi2 += res*res*w;
	  }
    }
  return chi2;

}

#include "fitsimage.h"

bool ImagePSF::FitPSF(PSFStarList &Stars)
{
  if (!refitPSF)
    {
      cout << "FitPSF : no psf refit requested for image " 
	   << reducedImage->Name() << endl;
      return true;
    }
  cout << " FitPSF: we have " << Stars.size() << " stars to fit the PSF " << endl;
  cout << " PSF half sizes = " << hSizeX << ' ' << hSizeY << endl;
  analyticPSF = ChooseAnalyticPSF(Cards().analyticKind);
  npar = analyticPSF->NPar();
  double seeing = reducedImage->Seeing();
  Vect startParams(npar); // initial uniform params
  analyticPSF->InitParamsFromSeeing(seeing,startParams);

  // in case some control output was requested:
  string outputDirectory; 
  if (Cards().outputDirectory != "") outputDirectory = Cards().outputDirectory;
  else outputDirectory = AddSlash(reducedImage->Dir());


  // load images.
  if(!FileExists(reducedImage->FitsWeightName())) {
    cout << "No weights. Exiting." << endl;
    return false;
  }
  FitsImage weight(reducedImage->FitsWeightName());
  // set weight of saturated pixels to 0.
  if (FileExists(reducedImage->FitsSaturName())) {
    FitsImage satur(reducedImage->FitsSaturName());
    weight *= (1-satur);
  }
  FitsImage image(reducedImage->FitsName());
  skylev = image.KeyVal("SKYLEV");

  // first fit of stars.
  for (PSFStarIterator i = Stars.begin(); i != Stars.end(); )
    {
      PSFStar &s = **i;
      s.SetPSFParams(startParams);
      if (s.FitPSFParams(image, weight, *this)) ++i;
      else 
	{
	  cout << " fit failure: erasing star at " 
	       << s.x << ' ' << s.y << ' ' << " fmax= " << s.Fluxmax() << endl;
	  i = Stars.erase(i);
	}
    }


  double mean,median,sigma, chi2;
  int ndof;
  Chi2Stat(Stars, mean, median, sigma, chi2, ndof, "independent profiles");

  // fit a smooth variation with position in frame of analytic PSF parameters, 
  double maxChi2 = median + Cards().nSigPSFChi2Cut * sigma;

  PSFStarList starsCopy(Stars); // just copy pointers, not objects; 

  if (!FitParametersVariation(starsCopy, maxChi2))
    {
      cout << " could not fit spatial variations " << endl;
      return false;
    }

  // Check that the found parametrization makes sense in the image corners
  bool paramsok = true;
  cout << " check parameters in the image corners" << endl;
  for (double x=0; x<nx+1; x+=nx)
  for (double y=0; y<ny+1; y+=ny)
    {
      cout << "parameters in (x,y) = (" << x << ',' << y << "):" << endl;
      Vect pars(PSFParams(x,y));
      pars.writeASCII(cout);
      if (!analyticPSF->CheckParams(pars))
	{
	  cout << " Unacceptable params at  (x,y) = (" << x 
	       << ',' << y << ")" << endl;
	  paramsok = false;
	}
    }

  if (!paramsok)
    {
      cout << " could not model a reasonable parameters spatial variation of the analytic PSF. stopping here " << endl;
      return false;
    }


  /* optimization loop of discretized residuals:
     1) fit all stars with current model
     2) discard outliers: if none and already have residuals and dchi2 small, break;
     3) compute residuals
     4) goto 1 (but no more than 5 iterations, because it takes time ...)
  */
     
  
  double oldChi2 = chi2;

  /* I first devised this routine as able to fit non linearity, but I
     now think it is useless. FitNonLinearity is the right routine.
     Since the code is there, I just disable it by setiing the degree
     to -1  */
  int nonLinDeg = -1;

  for (int grand_loop=0; grand_loop<2; ++grand_loop) // 2 iterations, with outlier pixel removal in between
    {
      int maxIterResidual[] = {5,3};

      for (int iterResiduals=0; iterResiduals < maxIterResidual[grand_loop]; ++iterResiduals)
	{

	  // debug
	  // double chi2_1 = ComputeChi2(Stars,image,weight);

	  NonLinModel *lastNonLin=0;
	  if (nonLinDeg >= 0)
	    lastNonLin = new NonLinModel(nonLinDeg, 0., 60000.);

	  // compute pixellized residuals.
	  
	  if (!StackResiduals(Stars, image, weight, lastNonLin))
	    {
	      cout << " could not fit pixellized residuals " << endl; 
	      return false;
	    }
	  
	  //debug
	  // double chi2_2 = ComputeChi2(Stars, image,weight);

	  if (Cards().allResidualsFileName != "")
	    for (unsigned k =0; k < residuals->NTerms(); ++k)
	      {
		char resFileName[64];
		sprintf(resFileName,"%s%s_%d_%d.fits",
			outputDirectory.c_str(),
			CutExtension(Cards().allResidualsFileName).c_str(),
			iterResiduals, k);
		residuals->Coeffs(k).writeFits(resFileName);
	      }
	  
	  if (nonLinDeg >= 0)
	    {
	      ApplyNonLinearities(image,*lastNonLin);
	      if (nonLinearity)
		{
		  (*nonLinearity) += (*lastNonLin);
		  delete lastNonLin;
		}
	      else nonLinearity = lastNonLin;
	    }

	  
	  cout << " refitting all stars with current PSF model " << endl;

	  int fitFailures = 0;

	  for (PSFStarIterator i = Stars.begin(); i != Stars.end();)
	    {
	      PSFStar &s = **i;
	      if (s.FitStarParams(image, weight, *this)) ++i;
	      else
		{
		  cout << " erasing star at " << s.x << ' ' << s.y 
		       << " for fit failure " << endl;
		  i = Stars.erase(i); 
		  fitFailures++;
		}
	    }
	  if (fitFailures) 
	    cout << " Erased " << fitFailures << " stars (fit failure) in the last fitting loop " 
		 << endl;

	  Chi2Stat(Stars, mean, median, sigma, chi2, ndof,"");
	  
	  maxChi2 = median + Cards().nSigPSFChi2Cut * sigma;
	  
	  // debug
	  // cout << " debug chi2 avant, milieu, fin " << chi2_1 << ' ' << chi2_2 << ' ' << chi2 << endl;

	  // discard outliers
	  int erasedCount = erase_outliers(Stars, maxChi2);

	  if (erasedCount)
	    {
	      cout << " Erased " << erasedCount << " stars (outliers) " << endl;
	      Chi2Stat(Stars, mean, median, sigma, chi2, ndof,"");// recompute stats
	    }

	  cout << " count , average chi2 ";
	  if (residuals) cout << "(global profile + discrete refinements) ";
	  else cout << "(global profile no discrete refinements) ";
	  cout << Stars.size() << ' ' << mean << ' ' << endl;
	  cout << " *** chi2/ndof " << chi2 << '/' << ndof << "=" << chi2/ndof << " **** dchi2 = " << oldChi2-chi2 << endl;

	  if (erasedCount + fitFailures == 0 && iterResiduals > 1 && fabs(oldChi2-chi2)<0.1) 
	    {
	      cout << " convergence on pixellized residuals reached " << endl;
	      break;
	    }
	  oldChi2 = chi2;
	}// end of loop without outlier pixel removal

      if (grand_loop == 0)
	{
	  int npix = RemoveOutlierPixels(Stars,image,weight,Cards().nSigPixResCut);
	  cout << " number of outlier pixels removed " << npix << endl;
	  if (npix ==0) break;
	}
    }

  // Done with the fit. outputs now:
  if (Cards().lastResidualsFileName != "")
    for (unsigned k =0; k < residuals->NTerms(); ++k)
      {
	char resFileName[64];
	sprintf(resFileName,"%s%s_%d.fits",
		outputDirectory.c_str(),
		CutExtension(Cards().lastResidualsFileName).c_str(),k);
	residuals->Coeffs(k).writeFits(resFileName);
      }
  

  // final output for grep in logs
  cout << " PSF_SUMMARY_name_nstars_chi2_ndof_chi2/dof " 
       << reducedImage->Name() << ' ' 
       << Stars.size() << ' '
       << chi2 << ' '
       << ndof << ' '
       << chi2/ndof << ' ' << endl;

  if (Cards().residualsTupleFileName != "")
    ResidualTuple(outputDirectory+Cards().residualsTupleFileName, Stars);
  if (Cards().residualsImageFileName != "")
    ResidualImage(outputDirectory+Cards().residualsImageFileName, Stars);

  Write(); // actually write the PSF
  /* .. and write the fitted stars (to be able to extract the fluxes 
     without refitting) */
  Gtransfo *wcs;
  if (!WCSFromHeader(image,wcs)) wcs = NULL;
  Stars.WriteTuple(reducedImage->Dir()+"psfstars.list",wcs, this);
  if (wcs) delete wcs;

  return true;
}

int ImagePSF::RemoveOutlierPixels(const PSFStarList &Stars,const Image &I, Image &W, const double NSigResCut) const
{
  int npix = 0;
  double chi2_cut = sq(NSigResCut);
  for (PSFStarCIterator it = Stars.begin(); it != Stars.end(); ++it)
    {
      const PSFStar &s = **it;
      double xc = s.PSFX();
      double yc = s.PSFY();
      double psfFlux = s.PSFFlux();
      int starti, endi, startj,endj;
      StampLimits(xc, yc, starti, endi, startj, endj);
      for (int j=startj; j <endj; ++j)
	for (int i=starti ; i < endi; ++i)
	  {
	    double w = W(i,j);
	    if (w<=0) continue;
	      
	    double psfVal = PSFValue(xc, yc, i, j);
	    double res = (I(i,j)-psfVal*psfFlux);
	    w = 1./(1./w+psfVal*psfFlux/gain);
	    if (res*res*w > chi2_cut)
	      {
		W(i,j) = 0;
		npix++;
	      }
	  }
    }
  return npix;
}


#include "sparsevect.h"
#include <ctime>


/* This is the key routine of this PSF modelling.
   It fits simultaneously:
      - the "PSF residuals" (pixels to be added to the analytic PSF)
      - the star fluxes
      - if requested the response non-linearity

   The reason to mix all these is that running separate minimizations 
   does not allow to reach rapidly the actual minimum, even if not fitting non-linearities.
   In the same spirit fitting non-linearities alone or together with fluxes does
   not converge.
   The code is reasonnably well tested and reasonnably fast, except in case of very poor
   image quality. To speed up in this specific case, we could perhaps fit a coarser
   pixellized model (since the PSF is "necessarily" oversampled) but it is a non-negligible amount of work.
   We may also give up spatial variations, because they become negligible when the seeing gets large.
*/

bool ImagePSF::StackResiduals(PSFStarList &Stars, 
			      const Image &I, const Image &W, NonLinModel *NonLin)
{
  int nstars = Stars.size();
  if (nstars == 0) return false;
  cout << " entering StackResiduals with half sizes = " << hSizeX << ' ' << hSizeY << ", nstars=" << nstars << endl;
  int hx = hSizeX;
  int hy = hSizeY;
  /* The resampling kernel size determines how larger we want the
   discretized PSF to be compared to the "user window" (the size used
   for fits) : 
    - if kernel size is 1 (nearest neighbor) increment by 0
    - if kernel size is 2, increment by 1,
    - if kernel size is 3, increment by 1
  This proved to be a bad idea, because the extra pixels,
  on the border are poorly defined, and cannot be used. 
  Poly2Image::PixValue assumes that these pixels are zero, as we do here.
*/
  //  int extraPix = ResampNCoeffs()/2;
  int extraPix = 0; // ResampNCoeffs()/2;
  int hxRes = hSizeX+extraPix;
  int nxRes = 2*hxRes+1;
  int hyRes = hSizeY+extraPix;
  int nyRes = 2*hyRes+1;

  Poly2Image tmpRes(0,0,nx,ny, Cards().spatialPixVarDeg, hxRes, hyRes);

  // image non-linearity stuff.
  unsigned nparNonLin = (NonLin)? NonLin->NPar() : 0;

  Vect hNonLin(nparNonLin);
  Vect hSkyNonLin(nparNonLin);

  if (NonLin)
  {
    FitsHeader head(reducedImage->FitsName());
    NonLin->ParamDerivatives(skylev, hSkyNonLin);
  }

  double maxPixVal = -1e30;

  int npix = tmpRes.NPix();
  int nterms = tmpRes.NTerms();
  int npix_coeffs = npix*nterms; 
  int startNonLin =  npix_coeffs+nstars;
  int constant_start = npix_coeffs+nstars + nparNonLin;
  unsigned ntot = constant_start+2*nterms; // why 2? 1 is for the spatially constant term.1 is for the constraint
  Mat a(ntot, ntot);
  Vect b(ntot);
  SparseVect h;

  Vect cx(ResampNCoeffs());
  Vect cy(ResampNCoeffs());
  Vect hpol(nterms);


  int ks = npix_coeffs;
  for (PSFStarCIterator it = Stars.begin(); it != Stars.end(); ++it, ++ks)
    {
      const PSFStar &s = **it;
      double xc = s.PSFX();
      double yc = s.PSFY();
      int ic = int(floor(xc+0.5));
      int jc = int(floor(yc+0.5));
      int minDi, maxDi, minDj, maxDj;
      ResampCoeffs(ic-xc,minDi,maxDi,cx);
      ResampCoeffs(jc-yc,minDj,maxDj,cy);
      tmpRes.Monomials(xc,yc, hpol);

      int starti, endi, startj,endj;
      starti = max(ic-hx,0);
      startj = max(jc-hy,0);
      endi = min(ic+hx+1,nx);
      endj = min(jc+hy+1,ny);
      

      double flux = s.PSFFlux();
      double totW = 0;


      for (int j=startj; j < endj; ++j)
	for (int i=starti; i < endi; ++i)
	  {
	    if (W(i,j) <= 0) continue;
	    h.zero();
	    for (int ip = 0; ip<nterms; ++ip) // loop on terms of the polynomial(s)
	    for (int jo = minDj; jo<= maxDj; ++jo)
	      {
		int jPix = j-jc+jo+hyRes;
		if (jPix <0 || jPix >= nyRes) continue;
		for (int io=minDi; io <= maxDi; ++io)
		  {
		    int iPix = i-ic+io+hxRes;
		    if (iPix<0 || iPix >=nxRes) continue;
		      {
			// position in the parameter vector of this "pixel"
			int hindex = npix*ip+iPix+jPix*nxRes; 
			h.set(hindex, hpol(ip)*cx(io-minDi)*cy(jo-minDj)*flux);			
		      }
		    }
		} // the pixels part of h is filled
	    double analyticPSFVal = 0;
	    double psfVal = PSFValue(xc, yc, i, j, NULL, NULL, &analyticPSFVal);
	    double expectedVal  = flux*psfVal;


	    /* derivative w.r.t non linearity parameters: these are
	       computed using expected value rather than observed one.
	       doing otherwise causes (very) large biases. See
	       README.ls_biases */
	    if (NonLin)
	      {
		NonLin->ParamDerivatives(expectedVal+skylev,hNonLin);
		hNonLin -= hSkyNonLin;
		for (unsigned ik=0; ik < nparNonLin; ++ik)
		  h.set(ik + startNonLin, -hNonLin(ik)); // minus sign is OK.
	      }

	    // derivative w.r.t flux (when flux is fitted)
	    h.set(ks, psfVal);

	    //derivative w.r.t the spatially constant term:
	    for (int ip = 0; ip<nterms; ++ip) // loop on terms of the residual polynomial(s)
	      {
		int hindex = constant_start+ip;
		h.set(hindex, hpol(ip)*flux*(1.-npix*analyticPSFVal));
	      }

	    h.sort(); // mandatory : h is filled in the wrong order to speed up computations

#if (false)
	    // CHECK for DEBUG
	    for (SVCIterator i1 = h.begin(); i1< h.end(); ++i1)
	      {
		 
		SVCIterator i2 = i1;
		++i2;
		if (i2 == h.end()) break;
		if (i1->index >= i2->index) abort();
	      }
#endif
		
	    /* We have to account for the contribution of the star to
               the photon noise, absent from the weight map. So we add it,
	       but using the model and *NOT* the data.See README.ls_biases.
	       Note that this introduces non-linearities.
            */
	    double pixWeight =  1./(1./ W(i,j)+expectedVal/gain);
	    if (pixWeight <= 0) continue;
	    totW += pixWeight;

	    double pixVal = I(i,j);
	    if (pixVal > maxPixVal) maxPixVal = pixVal;
	    double weightedRes =  pixWeight*(pixVal - expectedVal);

	    for (SVCIterator i1 = h.begin(); i1< h.end(); ++i1)
	      {
		int hindex1 = i1->index;
		double hval1 = i1->value;
		for (SVCIterator i2 = i1; i2< h.end(); ++i2)
		  {
		    a(hindex1, i2->index) += pixWeight*hval1*i2->value;
		  }
		b(hindex1) += weightedRes*hval1;
	      }
	  }// end loop on pixels for this star
      if (totW == 0) // this sometime happens, I don't quite get how, but it forbids to solve
	{
	  a(ks,ks) = 1;
	}
    } // end loop on stars

  /* Constraints: the sum of pixels is 0 at any place in the frame. means that the sum
  of all pixels of each "term" is zero 
  */
  for (int ip = 0; ip<nterms; ++ip)
    {
      int lambda = constant_start+nterms+ip;
      int ikmin = ip*npix;
      int ikmax = ikmin+npix;
      for (int ik = ikmin; ik < ikmax; ++ik)
	a(ik, lambda) = 1;
      b(lambda) = 0; // probably useless
    }



  /* check here if all diagonal terms of a are != 0 :we may
   miss a few pixels (sides or dead). If this happens the line and
   column also has to be 0 as well as the corresponding b
   component. replacing the diagonal null term of a by 1 should do a
   good job
  */
  for (int k = 0; k<npix_coeffs; ++k)
    if (a(k,k) == 0 ) 
      {
	cout << " ImagePSF::StackResiduals : unconstrained pix at index k = " << k << endl;
	a(k,k) = 1;
      }


  // the Lapack routines take hours if the matrix contains nan, so
  // check before proceeding
  for (unsigned k1=0; k1<ntot; ++k1)
    {
      for (unsigned k2=0; k2<ntot; ++k2)
	if (isnan(a(k1,k2)))
	  {
	    cout << " StackResiduals : the matrix contains nan's " << endl;
	    return false;
	  }
      if (isnan(b(k1)))
	{
	  cout << " StackResiduals : the vector contains nan's " << endl;
	  return false;
	}
    }

  //DEBUG check
  //  Mat cpa(a);
  // Vect cpb(b);

  //   b.writeASCII("b.dat");

  cout << " solving linear system (size = " << a.SizeX() << ")" << endl;
  //a.writeFits("a.fits");
  //b.writeASCII("b.dat");
  // since there are constraints, the matrix is NOT posdef.
  // So we cannot use the Cholesky factorization (dposv in lapack). But since it is
  // twice as fast as the next fastest factorization (dsysv in lapack), I cooked up a "block solver"
  // which uses Cholesky for the posdef part and hand made stuff for the remainder.

  clock_t tstart = clock();

  int info;
  if (NonLin == NULL)
    {
      info = cholesky_solve_quasi_posdef(a,b,constant_start,"U");
    }
  else
    {
      /* involved trick here: I wish to get the covariance matrix of
       the nonlin parameters.  rather than inverting the whole matrix,
       I solve the system for a few (=nparNonLin) RHS from which I can
       extract the small covariance matrix.  */
      Mat tempB(ntot, 1+nparNonLin);
      for (unsigned i=0; i<ntot; ++i) tempB(i,0) = b(i); // copy the "real" RHS

      // add the ones needed for inversion
      for (unsigned k=0; k <nparNonLin; ++k) tempB(startNonLin+k,1+k) = 1.; 
      // solve the system
      info = cholesky_solve_quasi_posdef(a,tempB,constant_start,"U"); 
      // place the result where the actual parameters are read downstream
      for (unsigned i=0; i<ntot; ++i) b(i) = tempB(i,0); 
      // collect non linearity parameters ...
      Vect tmp(nparNonLin);
      for (unsigned ik =0; ik < nparNonLin; ++ik) tmp(ik) = b(startNonLin+ik);
      NonLin->SetParams(tmp);
      // ... and  their covariance matrix
      Mat &nonLinCov = NonLin->Cov();
      for (unsigned i=0; i<nparNonLin; ++i)
	for (unsigned j=0; j<nparNonLin; ++j)
	  {
	    nonLinCov(i,j) = tempB(startNonLin+j, 1+i);
	  }
      NonLin->SetMaxPix(maxPixVal);
    } /* end of nonlin special processing. b contains the same value
	 as if not fitting non linearities. */

  if (info != 0)
      {
	cout << " could not solve for pixellized residuals " << endl;
	cout << " size of the matrix " << a.SizeX() << endl;
	cout << " npix = " << npix << "  nterms= " << nterms << endl;
	return false;
      }
      
  clock_t tend = clock();
  cout << " CPU for solving " << float(tend-tstart)/float(CLOCKS_PER_SEC) << endl;

  cout << " StackResiduals : max pix value =  " << maxPixVal << endl;

  // propagate fit results : propagate constants back into the psf pixels:
  for (int ip = 0; ip<nterms; ++ip)
    {
      double cst = b(constant_start+ip);
      int ikmin = ip*npix;
      int ikmax = ikmin+npix;
      for (int ik = ikmin; ik < ikmax; ++ik) b(ik) += cst;
    }

  // put the pixels into the psf pixels (Poly2Image class)
  tmpRes.SetCoeffs(b);
  Vect cst(nterms);

  if (!residuals) residuals = new Poly2Image(tmpRes);
  else *residuals += tmpRes;


  // propagate fit results :  update fluxes
  ks = npix_coeffs;
  for (PSFStarIterator it = Stars.begin(); it != Stars.end(); ++it, ++ks)
    {
      PSFStar &s = **it;
      s.SetPSFFlux(s.PSFFlux() + b(ks));
    }


  return true;
}


void ImagePSF::ApplyNonLinearities(Image &I, const NonLinModel& NonLin) const
{
  cout << " Applying non-linearity correction:" << endl;
  for (unsigned k=0; k < NonLin.NPar(); ++k) cout << NonLin.Params()(k) << ' ';
  cout << endl;
    
  double tfSky = NonLin.Value(skylev);
  Pixel *pend = I.end();
  for (Pixel *p = I.begin(); p < pend; ++p)
    {
      *p += NonLin.Value(*p+skylev)-tfSky;
    }
}


void ImagePSF::ResidualImage(const string &ResFitsName, 
			     const PSFStarList &Stars)const
{
  FitsImage image(reducedImage->FitsName());
  const FitsHeader &head = image;
  FitsImage res(ResFitsName, head);
  FitsImage stack(DirName(ResFitsName)+"/psf_res_stack.fits", 2*hSizeX+1, 2*hSizeY+1);
  for (PSFStarCIterator it = Stars.begin(); it != Stars.end(); ++it)
    {
      const PSFStar &s = **it;
      double x = s.PSFX();
      double y = s.PSFY();
      int starti, endi, startj,endj;
      StampLimits(x,y, starti, endi, startj, endj);
	for (int j=startj; j < endj; ++j)
	  for (int i=starti; i < endi; ++i)
	    {
	      double val  = PSFValue(x,y, i, j);
	      val *= s.PSFFlux();
	      res(i,j) = image(i,j) - val;
	      stack(i-starti, j - startj) += (image(i,j) - val)/s.PSFFlux();
	    }
    }
  stack *= (1./Stars.size());
}

#include <fstream>

void ImagePSF::ResidualTuple(const string &ResTupleName, 
			     const PSFStarList &Stars)const
{
  FitsImage image(reducedImage->FitsName());
  ofstream file(ResTupleName.c_str());
  file << "# xc : " << endl;
  file << "# yc : " << endl;
  file << "# flux : " << endl;
  file << "# ic : center pixel " << endl;
  file << "# jc : center pixel " << endl;
  file << "# obj: obj id" << endl;  
  file << "# i : offset w.r.t center " << endl;
  file << "# j : offset w.r.t center " << endl;
  file << "# fpsf :   " << endl;
  file << "# fimg : measured flux " << endl;
  file << "# end" << endl;

  int count  = 0;
  for (PSFStarCIterator it = Stars.begin(); it != Stars.end(); ++it)
    {
      const PSFStar &s = **it;
      double x = s.PSFX();
      double y = s.PSFY();
      int starti, endi, startj,endj;
      StampLimits(x,y, starti, endi, startj, endj);
      int icenter = int(floor(x+0.5));
      int jcenter = int(floor(y+0.5));
      count ++;
      for (int j=startj; j < endj; ++j)
	for (int i=starti; i < endi; ++i)
	  {
	    
	    double val  = PSFValue(x,y, i, j);
	    val *= s.PSFFlux();
	    //	      res(i,j) = image(i,j) - val;
	    file << x << ' ' 
		 << y << ' '
		 << s.PSFFlux() << ' '
		 << icenter << ' '
		 << jcenter << ' '
		 << count << ' ' 
		 << i - icenter << ' '
		 << j - jcenter << ' '
		   << val << ' ' 
		 << image(i,j) << ' '
		 << endl;

	  } // end loop on pixels
    }// end loop on objects
  file.close();
}

//#include "poly1withfixedvalue.h"

/*! This routine fits no-linearity of response of CCD.

It fits a polynomial non-linearity and refits fluxes at the same time.
The minimized chi^2 reads:

   sum_pixels   (F(p+sky) - F(sky) - (flux+dflux)*psf(i,j))**2 * w

   where: 
     - p is the (sky subtracted) pixel value
     -  F(x) = x + x**2*P(x)   where P is a polynomial
     - w is the pixel weight
          we use the weight map which only accounts for the sky and we
           add on the fly the contribution of the object light (from the model, **NOT** the data)
     -  One fundamental trick is to compute the derivatives of P at the **MODEL** (flux*psf) value
         and NOT at the observed pixel value. If you do not do that, the chi2 gradient is not zero
         on average when you set the parameters at their true value. The biases encountered if one
         fails to do that may be larger than a "reasonnable" non-linearity.

*/

bool ImagePSF::FitNonLinearity(const bool WriteResTuple) /* const */
{ 
  if (!analyticPSF)
    {
      cout << " no PSF stored .. giving up " << endl;
      return false;
    }

  int nonLinDeg = Cards().nonLinearityDegree;
  if (nonLinDeg < 0 )
    {
      cout << " fitting non linearity with degree = -1 : there is nothing to fit ! " << endl;
      return false;
    }

    // load images.
  if(!FileExists(reducedImage->FitsWeightName())) {
    cout << "No weights. Exiting." << endl;
    return false;
  }
  FitsImage weight(reducedImage->FitsWeightName());
  { // set weight of saturated pixels to 0.
    FitsImage satur(reducedImage->FitsSaturName());
    weight *= (1-satur);
  }
  FitsImage image(reducedImage->FitsName());
  skylev = image.KeyVal("SKYLEV");

  PSFStarList stars(AperSEStarList(reducedImage->StarCatalogName()));
  if (stars.size() == 0)
    {
      cout << " no stars .. giving up " << endl;
      return false;
    }
  int fitFailures = 0;
  

  cout << " refitting all stars " << endl;
  for (PSFStarIterator i = stars.begin(); i != stars.end();)
    {
      PSFStar &s = **i;
      if (s.FitStarParams(image, weight, *this)) ++i;
      else
	{
	  cout << " erasing star at " << s.x << ' ' << s.y 
	       << " for fit failure " << endl;
	  i = stars.erase(i); 
	  fitFailures++;
	}
    }
  if (fitFailures) 
    cout << " Erased " << fitFailures << " stars (fit failure) in the last fitting loop " 
	 << endl;


  // discard outlier stars
  int erasedCount;
  do
    {
        double mean,median,sigma, chi2_init;
	int ndof;
	Chi2Stat(stars, mean, median, sigma, chi2_init, ndof, "independent profiles");

	double maxChi2 = median + Cards().nSigPSFChi2Cut * sigma;  
	erasedCount = erase_outliers(stars, maxChi2);


	if (erasedCount)
	  {
	    cout << " Erased " << erasedCount << " stars (outliers) " << endl;
	    Chi2Stat(stars, mean, median, sigma, chi2_init, ndof,"");// recompute stats
	  }
    }
  while (erasedCount != 0);

  stars.sort(DecreasingFlux); // don't remember why .... 

  // clear outlier pixels.
  int npix = RemoveOutlierPixels(stars,image,weight,Cards().nSigPixResCut);
  cout << " Found and removed " << npix << " outlier pixels " << endl;
  

  NonLinModel nonLin(nonLinDeg, 0, 60000.);
  if (!StackResiduals(stars,image,weight, &nonLin)) return false;
  nonLin.Write((reducedImage->Dir()+"/nonlinpar.dat").c_str());

  // from here on, it is only diagnostic code..

  // write the "corrected" star list
  Gtransfo *wcs;
  if (!WCSFromHeader(image,wcs)) wcs = NULL;
  stars.WriteTuple((reducedImage->Dir()+"/psfstars2.list").c_str(),wcs,this);

  // write the "corrected" psf
  this->Write((reducedImage->Dir()+"/psf2.dat").c_str());


  // write the residual tuple if requested.
  ofstream restuple;


  if (WriteResTuple)
    {
      restuple.open((reducedImage->Dir()+"/restuple.list").c_str());
      
      restuple << "# xc :" << endl;
      restuple << "# yc : " << endl;
      restuple << "# i :" << endl;
      restuple << "# j :" << endl;
      restuple << "# w : weight " << endl;
      restuple << "# psfVal :" << endl;
      restuple << "# psfFlux : new psf flux" << endl;
      restuple << "# fmax :   peak star flux" << endl;
      restuple << "# imij :   old pix value" << endl;
      restuple << "# nimij :  new pix value" << endl;
      restuple << "# df :  new-old psf flux " << endl;
      restuple << "# nres : new residual" << endl;
      restuple << "# ores : old residual" << endl;
      restuple << "#end"<< endl;
    }


  for (unsigned k=0; k < nonLin.NPar(); ++k)
    cout << "param(" << k << ")= " << nonLin.Params()(k) << " +/- " 
	 << sqrt(nonLin.Cov(k,k)) << endl;
  
  double chi2New = 0;
  double chi2Old = 0;
  double sumDflux = 0;
  double sumFluxes = 0;
  for (PSFStarCIterator it = stars.begin(); it != stars.end(); ++it)
    {
      const PSFStar &s = **it;
      //      cout << " flux df/flux " << s.PSFFlux() << ' ' << b(ks)/s.PSFFlux() << endl;
      double xc = s.PSFX();
      double yc = s.PSFY();
      double psfFlux = s.PSFFlux();
      double oldPsfFlux = s.OldPSFFlux();
      double dflux = psfFlux - oldPsfFlux;
      sumDflux += dflux;
      sumFluxes += oldPsfFlux;
      

      int starti, endi, startj,endj;
      StampLimits(xc,yc,starti, endi, startj, endj);
      for (int j=startj; j < endj; ++j)
	for (int i=starti; i < endi; ++i)
	  {
	    double w = weight(i,j);
	    if (w<=0) continue;	    
	    double psfVal = PSFValue(xc,yc,i,j);
	    w = 1./(1./w+psfVal*psfFlux/gain);
	    double val1 = psfFlux*psfVal+skylev;
	    double pix= image(i,j) + (nonLin.Value(val1)-nonLin.Value(skylev));
	    double res = (pix- psfFlux*psfVal);
	    chi2New += res*res*w;
	    double ores = image(i,j)-oldPsfFlux*psfVal;
	    chi2Old += ores*ores*w;
	    if (WriteResTuple)
	      {
		restuple << xc << ' ' << yc << ' ' << i << ' ' << j << ' ' 
			 << w << ' ' << psfVal << ' ' << psfFlux << ' '
			 << s.Fluxmax() << ' ' << image(i,j) << ' '  << pix << ' '
			 << dflux << ' ' << res << ' ' << ores <<endl;
	      }
	  }
    }
  cout << " chi2 old new " << chi2Old << ' ' << chi2New << endl;
  cout << " sumfluxes sumDFlux " << sumFluxes << ' ' << sumDflux << endl;


  if (WriteResTuple) restuple.close();

  return true;
}




/* This MakePSF routine is in fact something like
ReducedImage::MakePSF. It was put here to avoid bindings between the
core poloka and the PSF code
*/
      
#include <sstream>
#include "nstarmatch.h"
#include "apersestar.h"
#include "wcsutils.h"
#include "fitsimage.h"
#include "reducedimage.h"
#include "listmatch.h"
#include "psfstar.h"
#include "astroutils.h" // IdentifyDeepField.

bool MakePSF(const string &ImageName, const bool RefitPSF,
	     const bool UseExternalCatalog)
{
  ReducedImage *ri = new ReducedImage(ImageName);
  ReducedImageRef reducedImage = ri;
  if (!RefitPSF && FileExists(PSFFileName(*ri))) 
    {
      cout << " MakePSF: image " << ri->Name() << " already has a psf file and no PSF refit requated. Nothing done " << endl;
      return true;
    }

  reducedImage->MakeAperCat();
  AperSEStarList starImageCat;

  if (UseExternalCatalog)
    {
      FitsHeader head(reducedImage->FitsName());
      string fieldName;
      if (IdentifyDeepField(head, fieldName))
	{
	  ostringstream starRefCatName;
	  starRefCatName << fieldName << ".list";
	  string wholeStarRefCatName = DbConfigFindCatalog(starRefCatName.str());
	  if (wholeStarRefCatName != "")
	    {
	      Gtransfo *readWcs;
	      TanPix2RaDec *wcs;
	      if (WCSFromHeader(head, readWcs) 
		  && ((wcs = dynamic_cast<TanPix2RaDec*>(readWcs))))
		{
		  BaseStarList starRefCat(wholeStarRefCatName);
		  cout << " read " << starRefCat.size() << " stars from " 
		       << wholeStarRefCatName << endl;

		  /* projection from sideral coordinates to tg plane
		     (with coordinates in degrees)
		  */
		  GtransfoLin id;
		  TanRaDec2Pix projection(id, wcs->TangentPoint());
		  Gtransfo *im2TgPlane = GtransfoCompose(&projection, wcs);

		  // project the star catalog to the tg plane:
		  starRefCat.ApplyTransfo(projection);
		  /* for information : non-destructive way of doing the 
		     same thing : 
		     TStarList projStarCatalog(starCat, projection);
		  */
		  AperSEStarList wholeImageCat(reducedImage->AperCatalogName());

		  StarMatchList *matches = ListMatchCollect(AperSE2Base(wholeImageCat), 
							    starRefCat,
							    im2TgPlane,
							    2./3600.);
		  cout << " found " << matches->size() << " matches " << endl;
		  delete im2TgPlane;
		  for (StarMatchIterator i=matches->begin(); i != matches->end(); ++i)
		    {
		      AperSEStar *s = dynamic_cast<AperSEStar*>(&(*(i->s1)));
		      if (s->Flag() || s->FlagBad()) continue;
		      starImageCat.push_back(s);
		    }
		  delete matches;
		}
	    }
	}
    }
  if (!UseExternalCatalog || starImageCat.size() == 0)
    {
      if (!reducedImage->MakeStarCat())
	{
	  cout << " ERROR : MakePSF: No PSF for " << reducedImage->Name() 
	       << "( missing star catalog) " << endl;
	  return false;
	}
      string catName = reducedImage->StarCatalogName();
      starImageCat.read(catName);
      cout << " MakePSF : using stars from " << catName << endl;
    }
  ImagePSF imagePSF(*reducedImage, RefitPSF);
  PSFStarList stars(starImageCat);
  return imagePSF.FitPSF(stars);
}
 
  
