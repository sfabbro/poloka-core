// This may look like C code, but it is really -*- C++ -*-
#ifndef KERNELFIT__H
#define KERNELFIT__H

#include <list>
#include <vector>


using namespace std;

#include "basestar.h"
#include "image.h"
#include "dimage.h"
#include "frame.h"
#include "countedref.h"

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

#include <vector>


struct XYPower  /* to handle and compute things like x**n*y**m */
{  
  /* was separated from other stuff because we need 2 of them : one for the spatial variations of the kernel,
   and one for the spatial variation of the background (see in OptParams) */
  int Degree;
  vector<int> Xdeg; // values of x exponant of monomials
  vector<int> Ydeg; // values of exponant for x and y of monomials
  size_t Nterms() const {return Xdeg.size();};

  /* the value of monomial of rank q (where q < Nterms) */
  double Value(const double X, const double Y, const size_t q) const;

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


class Mat;
class Vect;

/* The actual hanger for kernel fit data */

//! A class to fit a convolution kernel between 2 images by least squares.
class KernelFit  :  public RefCount 
{
  public :
//! pointer to 'best' image (smaller seeing).
  const Image *BestImage; //!
//! pointer to 'worst' image (larger seeing).
  const Image *WorstImage; //!
  //!
  const Image *BestImageWeight;
  const Image *WorstImageWeight;

//! data frame
  Frame DataFrame; //!
//! the value of the sky of the best image.
  double BestImageBack;
//! the value of the sky of the worst image.
  double WorstImageBack;
//! the value of the sky variance of WorstImage.
  double SkyVarianceWorstImage;
  //! the gain of WorstImage
  double WorstImageGain;

//! the kernel integral at the image center (photometric ratio)
  double KernAtCenterSum; 
  vector<Kernel> Kernels; //! /* the base used to build the 'best' kernel */
//! the list of stamps in the 'best' image

  StampList  *BestImageStamps;

  //! a pointer to the 'best' image convolved with the solution
  Image *ConvolvedBest; //!
  DImage* convolutions;  //!
  DImage* backStamps;  //!
  OptParams optParams;
  //! a pointer to the differential background image
  Image *WorstDiffBkgrdSubtracted;


  int HKernelSizeX() const { return optParams.HKernelSize;}
  int HKernelSizeY() const { return optParams.HKernelSize;}

  size_t mSize; /* the size of the matrix m */
#if 0
  Mat m;  /* the least-squares matrix (the second derivative of chi2 w.r.t fitted parameters */ // we don't save it  
  Vect b;  /* the normal equations (i.e. grad(chi2) = 0) right hand side */ // we don't save it 
#endif

  vector<double> solution; //[mSize]/* the weights of various kernels */
  vector<double> diffbackground; //! /* the differential background coefficient when fitted separately */

/* routines to compute an index in the solution array ... internal cooking */

  int KernIndex(int iKernel, int iSpatial) const { return optParams.KernVar.Nterms()*iKernel + iSpatial;};
  int BackIndex(int ib) const { return ib+Kernels.size()*optParams.KernVar.Nterms();}


/* Fit the differential background separately from the kernel */
  int FitDifferentialBackground(const double NSig);

/* fitted differential background value */
  double BackValue(const double&x, const double &y) const ; 

  void DeleteConvolvedBest() { if (ConvolvedBest) delete ConvolvedBest; ConvolvedBest = NULL;}
  void DeleteWorstDiffBkgrdSubtracted() { if (WorstDiffBkgrdSubtracted) delete WorstDiffBkgrdSubtracted; WorstDiffBkgrdSubtracted = NULL;}

  KernelFit() { 
    BestImage = NULL;
    WorstImage = NULL; 
    BestImageWeight = NULL;
    WorstImageWeight = NULL; 
    BestImageStamps = NULL; 
    ConvolvedBest = NULL; 
    WorstDiffBkgrdSubtracted = NULL;
    convolutions = NULL;
    BestImageBack = 0;
    WorstImageBack = 0;
    nstamps=0;
  }
  
  ~KernelFit() {
    if (convolutions) delete [] convolutions;
    if (BestImageStamps) delete BestImageStamps;
    if (ConvolvedBest) delete ConvolvedBest;
    if (WorstDiffBkgrdSubtracted) delete WorstDiffBkgrdSubtracted;
  }

  
/* fitting routines */
void KernelsFill(); /* compute the kernels constituting the basis */
void AllocateConvolvedStamps(); /* trick to avoid allocation and deallocation of a big memory chunk. 
                                 triggered when needed */
void DeallocateConvolvedStamps();

//! computes least square matrix and vector for one stamp
  //void OneStampMAndB(const Stamp &Astamp, vector<double> &StampM, vector<double> &StampB);
  void OneStampMAndB(const Stamp &Astamp, Mat &StampM, Vect &StampB);

//! actually solves the system by calling linear algebra efficient routines
  bool Solve(const Mat &M, const Vect &B);

  double StampChi2(Stamp &stamp, double VSky, double InvGain); /* involves a kernel computation and a stamp convolution */
  double chi2;
  int NStampsUsed() const {return nstamps;}

/* computes chi2 contributions of stamps and applies median filtering for outliers removal. iterates until
the number of stamps involved stabilizes. */
  void FilterStamps(Mat &M, Vect &B); 
  //! final wrapper that calls various routines to fill matrix and vector, and then solve.
int DoTheFit(const BaseStarList &List, double &BestSeeing, double &WorstSeeing);
int DoIt(const BaseStarList &List, double &BestSeeing, double &WorstSeeing);

/* computes the kernel at the given location according to the present solution */
void KernCompute(Kernel &Result, const double X, const double Y) const; 


void BestImageConvolve(int UpdateKernStep = 100);

Image *VarianceConvolve(const Image &Source, int UpdateKern = 100);

  
  //! read object contents from a stream
  void read(std::istream& stream);
  
  //! read object contents to a stream
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
  
#ifndef SWIG
  //  friend ostream& operator << (ostream &stream, const KernelFit &s)
  //  { s.dump(stream); return stream;}
  friend ostream& operator << (ostream &stream, const KernelFit &s)
  { s.write(stream); return stream;}
  friend istream& operator >> (istream &stream, KernelFit &s)
  { s.read(stream); return stream;}
#endif

private:
  int nstamps;

};


typedef CountedRef<KernelFit> KernelFitRef;


#endif /* KERNELFIT__H */
