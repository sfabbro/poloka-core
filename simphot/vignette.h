#ifndef VIGNETTE__H
#define VIGNETTE__H


#include <string>

#include "reducedimage.h"
#include "gtransfo.h"
#include "countedref.h"


#include "intframe.h"
#include "pixelblock.h"
#include "intpoint.h"

class SimPhotFit;
class Array4D;
class Mat;
class Vect;
class ImagePSF;




class Vignette : public RefCount {
 private :
  ReducedImageRef ri;
  double mmjd, mjd, seeing, expTime;
  bool couldFitFlux;
  SimPhotFit &simPhotFit;
  GtransfoRef vignette2Model; // image to model transfo (stamp coordinates)
  GtransfoRef model2Vignette; // image to model transfo (stamp coordinates)

  Point posInStamp; // coordinates of the SN in the local coordinates
  IntPoint intPos;  // nearest pixel

  
  double photomRatio; // definition : ref*photomRatio = this
  double flux,sky;
  double chi2;
  int nterms;

  const ImagePSF *imagePSF;

  // pixels handling 
  /* may be the IntFrames could be replaced by the corresponding stuff in PixelBlock's */
  IntFrame convolvedStampLimits;
  IntFrame stampLimits;
  PixelBlock imagePix, weightPix, residuals;

  //
  PixelBlock kernel;


 public:

  Vignette(SimPhotFit &SPF,const ReducedImage &Current);

  string Name() const { return ri->Name();}

  bool SetGeomTransfos();
  bool SetKernel();
  int HalfKernelSize() const { return (kernel.xmax-kernel.xmin)/2;}
  int CanDo() const;

  void ComputeModelLimits(Frame &ModelFrame) const;

  void FillAAndB(Mat &A, Vect &B, const int ToDo);

  void UpdateResiduals();

  //! assumes that residuals are up to date. returns number of "killed" pixels
  int KillOutliers(const double NSigCut);

  double Chi2() const { return chi2;}
  int NTerms() const { return nterms;}


  void Write(const string &Directory) const;

  double GetSky() const {return sky;}
  double GetFlux() const {return flux;}
  void SetSky(const double &Val) { sky = Val;}
  void SetFlux(const double &Val) { flux = Val;}


  double MJD() const { return mjd;}
  double MMJD() const { return mmjd;}
  double Seeing() const {return seeing;}
  double ExpTime() const { return expTime;}

 private :

  void ComputeGalaxyDerivatives(Array4D &);

  bool ReadPixels();

  bool GetPSF(PixelBlock &PSFPixels, PixelBlock *PSFDerX = NULL,
	      PixelBlock *PSFDerY=NULL) const;


  typedef CountedRef<Vignette> VignetteRef;


  friend bool IncreasingDate(const VignetteRef &V1, const VignetteRef &V2);
};


typedef CountedRef<Vignette> VignetteRef;
//class VignetteRef : public CountedRef<Vignette> {};

inline bool IncreasingDate(const VignetteRef &V1, const VignetteRef &V2)
{
  return (V1->MMJD() < V2->MMJD());
}




#include <list>




typedef list<VignetteRef> VignetteList;
typedef VignetteList::const_iterator VignetteCIterator;
typedef  VignetteList::iterator VignetteIterator;



class Array4D;

void ConvolveModel(const PixelBlock &M, const Array4D &Coeffs, 
		   PixelBlock &Result);




#endif /* VIGNETTE__H */
