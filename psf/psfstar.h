#ifndef PSFSTAR__H
#define PSFSTAR__H



#include "apersestar.h"
#include "matvect.h"

class Image;
class ImagePSF;
class Mat;
class Vect;

class PSFStar : public AperSEStar {

 protected :
  double psfX, psfY, psfFlux;
  double psfChi2;
  double oldPsfFlux;
  Vect psfParams;
  Mat  psfParamsWeight;
  Mat  xyfCov;
  
 public :

  PSFStar(const AperSEStar &A);
  PSFStar(const double X, const double Y, const double Flux);
  bool FitPSFParams(const Image &I, const Image &W, const ImagePSF &PSF);
  bool FitStarParams(const Image &I, const Image &W, const ImagePSF &PSF);

  void SetPSFParams(const Vect &PSFParams);
  const Vect& PSFParams() const {return psfParams;}
  const Mat& PSFParamsWeight() const {return psfParamsWeight;}
  const Mat& XYFCov() const{ return xyfCov;}
  double PSFChi2() const { return psfChi2;}
  void SetPSFChi2(const double &Chi2) {psfChi2 = Chi2;}
  double PSFX() const { return psfX;}
  double PSFY() const { return psfY;}
  double PSFFlux() const { return psfFlux;}

  double OldPSFFlux() const { return oldPsfFlux;}

  void SetPSFFlux(const double &Val) { psfFlux = Val;}



};


class Gtransfo;
#include "starlist.h"

class PSFStarList : public StarList<PSFStar> {

 public :
  PSFStarList () {} ;
  PSFStarList(const AperSEStarList &L);

  void WriteTuple(const string &FileName, const Gtransfo *Wcs= NULL, const ImagePSF* PSF=NULL) const;

  bool ReadTuple(const string &FileName);


};

typedef PSFStarList::const_iterator PSFStarCIterator;
typedef PSFStarList::iterator PSFStarIterator;



#endif /* PSFSTAR__H */
