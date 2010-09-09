// This may look like C code, but it is really -*- C++ -*-
#ifndef KERNELFIT__H
#define KERNELFIT__H

#include <list>
#include <vector>


using namespace std;

#include "dimage.h" // needed for vector<Kernel>

class StarMatchList;
class ImagePair;
class Image;

/*! \class KernelFit
    \brief Kernel fitting by least squares. 

Implementation of Ch. Alard's method for fitting (slowly variable) 
convolution kernels. see astro-ph 9903111 & astro-ph 9712287 
(where most of the equations of page 2 are wrong).
 Main differences with the "isis" code provided by Alard himself:
\arg if the integral of the kernel is requested to be constant over the image,
        we use Lagrange multipliers instead of something I could not understand.
        (paper from 97 page 3 first column).
\arg for the detection of outliers, we compute the chi2 of all stamps, 
     apply median filtering, refit, and iterate until stabilization
      of the number of stamps. This is in fact fairly fast because we have routines to 
      add or remove one stamp from the fit.
*/        


//  Short summary of classes involved:
//   The hanger for kernel fitting data is called KernelFit.
//  it contains among dozens of items a vector of Kernels  (vector<Kernel>).
// Since Kernel's are in double precision, they are NOT derived from Image's
// but from a new type DImage (it seemed to render things rather complicated
// to turn Image into a template class (on the pixel type used for calculations)).
// So in the end we have :
#ifdef JUST_FOR_DOCUMENTATION

//   class KernelFit { 
//   vector<Kernel> Kernels; 
//   BestImageStampList *StampList;
//   ...
//   };

#endif



struct XYPower  /* to handle and compute things like x**n*y**m */
{  
  /* was separated from other stuff because we need 2 of them : one for the spatial variations of the kernel,
   and one for the spatial variation of the background (see in OptParams) */
  int Degree;
  vector<int> Xdeg; // values of x exponant of monomials
  vector<int> Ydeg; // values of exponant for x and y of monomials
  unsigned nterms;
  unsigned Nterms() const {return nterms;};

  /* the value of monomial of rank q (where q < Nterms) */
  double Value(const double X, const double Y, const unsigned q) const;

  //!
  void ApplyToImage(Image &I, double Factor, const vector<double> &ParamVal) const;

  /* default constructor: Value(x,y,0) will return 1. */
  XYPower() { SetDegree(0);};
  XYPower(int Degree) { SetDegree(Degree);};
  void SetDegree(const int Degree);
  ~XYPower() {};
  
  //! read the polynomial from a stream
  void read(std::istream& stream);
  
  //! write the polynomial to a stream
  void write(std::ostream& stream) const;
  
  void dump(std::ostream &stream = std::cout) const{
    stream << "XYPower Xdeg ";
    for(size_t i=0;i<Xdeg.size() ;++i)
      stream << Xdeg[i] << " ";
     stream << "Ydeg ";
    for(size_t i=0;i<Ydeg.size() ;++i)
       stream << Ydeg[i] << " ";    
  }
#ifndef SWIG
  friend ostream& operator << (ostream &stream, const XYPower &s)
  { s.write(stream); return stream;}
  friend istream& operator >> (istream &stream, XYPower &s)
  { s.read(stream); return stream;}
#endif
};

class OptParams {
public :
  int HKernelSize; /*  actual size = 2*HKernelSize+1 */

  /* kernel basis */
  int NGauss;
  vector<double> Sigmas;// [NGauss]
  vector<int> Degrees; // [NGauss]
  double NSig; // size (in pixel) of kernel in number of sigmas of guessed size.
  
  /* Kernel variations */
  XYPower KernVar;

  /* background */
  XYPower BackVar;

  /* background separately fitted*/
  XYPower SepBackVar;

  int HStampSize;

  int MaxStamps;
  bool UniformPhotomRatio;
  bool OrthogonalBasis;
  bool SubtractNoise;

  OptParams(); /* default values */
  void OptimizeSizes(double BestSeeing, double WorstSeeing);
  void OptimizeSpatialVariations(const int NumberOfStars);
  int StampSize() const { return 2*HStampSize+1;}
  /* DOCF returns the size of the convolved stamps */
  int   ConvolvedSize() const { return 2*(HStampSize - HKernelSize) + 1;}
  
  //! read object contents from a stream
  void read(std::istream& stream);
  
  //! write object contents to a stream
  void write(std::ostream& stream) const;
  
  void dump(std::ostream &stream = std::cout) const {
    stream << "OptParams " << endl;
    stream << " KernVar " <<  KernVar << endl;
    stream << " BackVar " <<  BackVar << endl;
    stream << " SepBackVar " <<  SepBackVar << endl;
  }
#ifndef SWIG
  //  friend ostream& operator << (ostream &stream, const OptParams &s)
  //  { s.dump(stream); return stream;}
  friend ostream& operator << (ostream &stream, const OptParams &s)
  { s.write(stream); return stream;}
  friend istream& operator >> (istream &stream, OptParams &s)
  { s.read(stream); return stream;}
#endif 
};

class Stamp;
class StampList;
class Mat;
class Vect;

/* The actual hanger for kernel fit data */

//! A class to fit a convolution kernel between 2 images by least squares.
class KernelFit 
{

  vector<Kernel> Basis; //! /* the base used to build the 'best' kernel */
//! the list of stamps in the 'best' image

  OptParams optParams;
  vector<double> diffbackground; //! /* the differential background coefficient when fitted separately */

  unsigned mSize; /* the size of the matrix m */
  int nstamps;
  double kernAtCenterSum;
  bool fitDone;

  friend class FitWorkSpace;

public :
  vector<double> solution; // the coefficients of the linear combination

  public :

  KernelFit() {mSize=0; nstamps=0; kernAtCenterSum=0; fitDone = false;  }

  int HKernelSizeX() const { return optParams.HKernelSize;}
  int HKernelSizeY() const { return optParams.HKernelSize;}

  //! is it useful ? 
  const OptParams & KernelOptParams() const { return optParams;}

  int NStampsUsed() const {return nstamps;}

  unsigned NParams() const { return mSize;}

  bool FitDone() const { return fitDone;}

  //! Photometric Ratio
  double KernAtCenterSum() const { return kernAtCenterSum;}

  //! the same
  //! Photometric Ratio
  double PhotomRatio() const { return kernAtCenterSum;}


  double Chi2() const { return chi2;}

  //! final wrapper that calls various routines to fill matrix and vector, and then solve.
  int DoTheFit(ImagePair &ImPair);

  /*! Computes the kernel at the given location according to the present solution */
  void KernCompute(Kernel &Result, const double X, const double Y) const; 

  /*! Same as above, but allocates it beforehand */
  void KernAllocateAndCompute(Kernel &Result, const double X, const double Y) const; 


  void ImageConvolve(const Image &In, Image &Out,int UpdateKernStep = 100);

  void VarianceConvolve(const Image &Source, Image &Out, int UpdateKern = 100);

  //!
  void read(const std::string &FileName);
  
  //! read object contents from a stream
  void read(std::istream& stream);

  //! write object contents to a file
  void write(const std::string& FileName) const;

  //! write object contents to a stream
  void write(std::ostream& stream) const;
  
  void dump(std::ostream &stream = std::cout) const{
    stream << "KernelFit ";
    stream << optParams << std::endl;
    stream << "solution ";
    for(unsigned int i=0;i<solution.size() ;++i)
      stream << solution[i] << " ";
    stream << endl;
    stream << "diffbackground ";
    for(unsigned int i=0;i<diffbackground.size() ;++i)
      stream << diffbackground[i] << " ";
    stream << endl;
  }



private :

/* routines to compute an index in the solution array ... internal cooking */

  int KernIndex(int iKernel, int iSpatial) const { return optParams.KernVar.Nterms()*iKernel + iSpatial;};
  int BackIndex(int ib) const { return ib+Basis.size()*optParams.KernVar.Nterms();}


/* Fit the differential background separately from the kernel */
  int FitDifferentialBackground(ImagePair &ImPair, const double NSig);

/* simultaneously fitted differential background value */
  double BackValue(const double&x, const double &y) const; 
/* separatly fitted differential background value */
  double SepBackValue(const double&x, const double &y) const; 
  
/* fitting routines */
  void BasisFill(Mat *BasisTransfo); /* compute the kernels constituting the basis */


  //! computes least square matrix and vector for one stamp
  //void OneStampMAndB(const Stamp &Astamp, vector<double> &StampM, vector<double> &StampB);
  void OneStampMAndB(const Stamp &Astamp, const Image &WorstImage, Mat &StampM, Vect &StampB);

  void OneStampNoiseMatrix(const Stamp &AStamp, const Image &WorstImage, Mat &M) const;

//! actually solves the system by calling linear algebra efficient routines. Center provided for printing purposes
  bool Solve(const Mat &M, const Vect &B, const Point &Center);

  double StampChi2(Stamp &stamp, const Image &WorstImage); /* involves a kernel computation and a stamp convolution */
  double chi2;

  //! for fit studies
  void ParameterGroups(Mat &Groups) const;



public :

#ifndef SWIG
  //  friend ostream& operator << (ostream &stream, const KernelFit &s)
  //  { s.dump(stream); return stream;}
  friend ostream& operator << (ostream &stream, const KernelFit &s)
  { s.write(stream); return stream;}
  friend istream& operator >> (istream &stream, KernelFit &s)
  { s.read(stream); return stream;}
#endif

  
};



#endif /* KERNELFIT__H */
