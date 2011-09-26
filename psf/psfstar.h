#ifndef PSFSTAR__H
#define PSFSTAR__H



#include "apersestar.h"
#include "matvect.h"

class Image;
class ImagePSF;
class Mat;
class Vect;


class PSFStar : public BaseStar {

 private :

  bool do_the_fit(const Image &I, const Image &W, 
		  const ImagePSF &PSF, const bool FitPos, const bool FitParams);



 protected :
  double eflux;
  double fluxmax;
  double psfX, psfY;
  double psfChi2;
  Vect psfParams;
  Mat  psfParamsWeight;
  Mat  fxyCov;
  
 public :
  double oldPsfFlux;

 public :

  PSFStar();

  PSFStar(const SEStar &A);
  PSFStar(const double X, const double Y, const double Flux);

  double Fluxmax() const { return fluxmax;}
  double& Fluxmax() { return fluxmax;}

  bool FitAllParams(const Image &I, const Image &W, 
		    const ImagePSF &PSF);

  bool FitStarParams(const Image &I, const Image &W, 
		     const ImagePSF &PSF);


  void SetPSFParams(const Vect &PSFParams);
  const Vect& PSFParams() const {return psfParams;}
  const Mat& PSFParamsWeight() const {return psfParamsWeight;}
  const Mat& FXYCov() const{ return fxyCov;}
  double PSFChi2() const { return psfChi2;}
  void SetPSFChi2(const double &Chi2) {psfChi2 = Chi2;}
  double PSFX() const { return psfX;}
  double PSFY() const { return psfY;}
  double EFlux() const { return eflux;}


  //I/O's
  void writen(ostream &s) const;
  std::string WriteHeader_(std::ostream &s, const char* i) const;

  static PSFStar* read(fastifstream & Rd, const char *Format); 
  void read_it(fastifstream & Rd, const char *Format); 


};


class Gtransfo;
#include "starlist.h"

class PSFStarList : public StarList<PSFStar> {

 public :
  PSFStarList () {} ;
  PSFStarList(const std::string &FileName) { read(FileName);}

  PSFStarList(const AperSEStarList &L);

  void WriteTuple(const string &FileName, const Gtransfo *Wcs= NULL, const ImagePSF* PSF=NULL) const;

  bool ReadTuple(const string &FileName);


};

typedef PSFStarList::const_iterator PSFStarCIterator;
typedef PSFStarList::iterator PSFStarIterator;



#endif /* PSFSTAR__H */
