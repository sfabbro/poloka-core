#ifndef SIMPHOTFIT__H
#define SIMPHOTFIT__H

#include <map>

#define SIMPHOTFIT_CHECK_BOUNDS


#include "reducedimage.h"
#include "vignette.h"
#include "lightcurvefile.h"
#include "model.h"
#include "dimage.h" // for Kernel class
#include "matvect.h"

#define FIT_GALAXY 1
#define FIT_SKY 2
#define FIT_FLUX 4
#define FIT_POS 8

#define FIT_ALL 15

using namespace std;


class CalibratedStar;

class SimPhotFit : public Model
{
 private:
  ObjectToFit objToFit; // could perhaps be a reference...
  const LightCurveFile &lcFile;
  VignetteList vignetteList;
  int maxKernelSize;
  int toDo;
  Mat A; // contains either the weight param matrix, or towards the end, the covariance, and sometimes even the Cholesky-factorized weight matrix
  Vect B;
  Mat NightMat;      // see fillNightMat
  Mat A_with_posfitted;  // contains  the weight param matrix of the fit without fixed position


 public :
  SimPhotFit(const ObjectToFit &O, const LightCurveFile &LCFile);

  const ObjectToFit &ObjToFit() const { return objToFit;}

  bool DoTheFit();
  bool OneMinimization(const int ToDo, const int MaxIter, 
		       const double DeltaChi2);

  bool BuildVignettes();

  void FindModelBoundaries();
  
  bool Write(const string &Directory, const bool WriteVignettes, const bool WriteMatrices);

  void WriteTupleHeader(ostream &Stream, const int NStars) const;
  
  void WriteTupleEntries(ostream &Stream, const CalibratedStar &CStar) const;

  void fillNightMat();


 private :


  void CumulateChi2(double &Chi2, int &NDof, const bool Print=false) const;
  void AssignIndicesAndToDo(const int CurrentToDo);
  void UpdateResiduals();
  bool OneIteration(const int CurrentToDo);
  void DispatchOffsets(const Vect& Offsets, const double Fact=1., 
		       const bool Verbose=true);

  
  // what concerns the handling of fitted parameters

 private:

  // cannot be "const Vignette*" because we use these map to update sky and flux
  typedef map<const Vignette*,int> VignetteMap;
  VignetteMap skyMap;
  VignetteMap fluxMap;
  VignetteMap toDoMap;

  int posIndex;
  int matSize;
  int nParamGal;

 public:
  int SkyIndex(const Vignette*) const;
  int FluxIndex(const Vignette *) const;
  int PosIndex() const {return posIndex;}
  
};

#include "polokaexception.h"

class SimPhotFitException : public PolokaException
{

 public :
  SimPhotFitException(const string &msg) : PolokaException(msg) {}

};


int DecodeFitType(const int FitType);

/* This routine computes the PSF and position derivatives.
   If PSFPixels has a non zero size, then it is the "frame"
over which the PSF is calculated. If PSFPixels has a zero size,
the area over which the PSF was computed is used.
*/
class ImagePSF;


#endif /* SIMPHOTFIT__H */


