#ifndef DETECTION__H
#define DETECTION__H

#include <string>
#include "basestar.h"

class Detection : public BaseStar
{  
  CLASS_VERSION(Detection,1);
  #define Detection__is__persistent 
  
 private:
  double eFlux;
  double sig2Noise;
  double aperFlux;
  double sig2NoiseCv;
  double xs; // start position
  double ys; // start position
  double varXX; // pos error. ordering xx,yy,xy relevant...
  double varYY; // pos error
  double varXY; // pos error
  int area; // number of pixels  used in aperflux
  int nBad; // number of bad pixels (i.e. w=0)
  double distBad; // distance to nearest bad pixel
  double mxx; // shape param
  double myy;// shape param
  double mxy;// shape param
  // Reference related scores
  double ra;
  double dec;
  double vRaRa, vDecDec , vRaDec; // ordering important.
  double fluxRef; // flux under the SN
  double prctInc; // flux/fluxref if fluxref>0 100. if not
  double xObj, yObj; // coordinates of object on ref
  double fluxObjRef; // flux of nearest object
  double distObjRef; // distance to the latter  
  //  double fwhmRef; // FWHM of nearesest object (obsolete)
  // double shapeRef; // flux/fluxmax of nearest object on ref (obsolete)
  double localback ; // local background as computed by detector
 
    //! temporary, for IO's
  void read_it(istream& r, const char* Format);


  public:
  Detection(const double X=0, const double Y=0, const double Flux=0);

  double LocalBack() const {return localback;}
  double &LocalBack() { return localback;}


  double EFlux() const {return eFlux;}
  double &EFlux() { return eFlux;}

  double AperFlux() const {return aperFlux;}
  double &AperFlux() {return aperFlux;}

  double Mxx() const { return mxx;}
  double &Mxx() { return mxx;}

  double Mxy() const { return mxy;}
  double &Mxy() { return mxy;}
  
  double Myy() const { return myy;}
  double &Myy() { return myy;}

  double VarXX() const { return varXX;}
  double &VarXX() { return varXX;}

  double VarXY() const { return varXY;}
  double &VarXY() { return varXY;}
  
  double VarYY() const { return varYY;}
  double &VarYY() { return varYY;}

  double VarRaRa() const { return vRaRa;}
  double &VarRaRa() { return vRaRa;}

  double VarRaDec() const { return vRaDec;}
  double &VarRaDec() { return vRaDec;}
  
  double VarDecDec() const { return vDecDec;}
  double &VarDecDec() { return vDecDec;}

  int Area() const {return area;}
  int& Area() { return area;}

  int NBad() const { return nBad;}
  int &NBad() { return nBad;}

  double Sig2Noise() const { return sig2Noise;}

  double Ra() const {return ra;}
  double &Ra() { return ra;}

  double Dec() const {return dec;}
  double &Dec() { return dec;}

  friend class DetectionProcess;

  // IO's
  std::string WriteHeader_(ostream & stream, const char*i) const;
  static Detection* read(istream& r, const char* Format);
  void writen(ostream &s) const ;
};


#include "starlist.h"


typedef StarList<Detection> DetectionList;
typedef DetectionList::iterator DetectionIterator;
typedef DetectionList::const_iterator DetectionCIterator;


typedef CountedRef<Detection> DetectionRef;


BaseStarList *Detection2Base(DetectionList *D);


class Image;
class ReducedImage;

class DetectionProcess { 
  private :
    double seeing; // actual image seeing
  double sigFilter; // sigma of gaussian filetr used for convolution.
    double cvImageNormalizedSig; // sigma of convolver image normalized to convolved variance. should be 1 in the absence of pixel correlations.
// average of convolver image normalized to convolved variance. should be 0.
    double cvImageNormalizedMean; 

     double radBackMin;// radius of inner disk avoided in the  area used to estimate background
    double radBackMax;// half side of square area used to estimate background

    Image *image; // the image
    Image *imagecv; // ... convolved with a gaussian filter
    Image *weight; // weight image
    string imageName;
    string weightName;
    // detection  cuts : 
    double nSigDet; // cut on convolved image
    double nSig; // final cut on flux
    double radAper;
    double radWeight;
    double radBad;
    double minDistDouble;
    

 public:
      //!
      DetectionProcess(const string &ImageName, 
				   const string &WeightName, 
		       const double Seeing, const double SigFilter);

    //! run the default detection and apply cuts
    void DoDetection(DetectionList&);

    //! set scores from positions (e.g. to get detections scores if below cuts)
    void DetectionScoresFromPositions(const BaseStarList &In, 
				      DetectionList &Out, 
				      bool FixedPositions = false);

    //! Detections  have a few ref related scores. This routine sets them ...
    /*! ...but DetectionProcess should be created  with the ref images names */
    void SetScoresFromRef(DetectionList &List, const ReducedImage &Ref);

    ~DetectionProcess();

 private:

    void FluxFromPos(const Point &Where, double &Flux);

    void SetCuts(const int ToPrint = 0);

    void ImagesConvolve();

    void SaveNormalizedSig() const;
    bool GetNormalizedSig();

    //! detections above threshold
    void FirstDetections(DetectionList &) const;
    void DetectionPosition(Detection &Det) const;
    void SetDetectionScores(Detection &Det) const;
    void RefineDetectionPosition(Detection &Det) const;
    double LocalImageBackground(const int i, const int j) const;



};

#include <vector>

  
class MatchedDetection : public Detection
{
 private:
  vector<DetectionRef> others;

 public:
  MatchedDetection() {}
  MatchedDetection(const Detection &D) : Detection(D) {};

  void AssociateDet(Detection *OtherDet )
    {others.push_back(OtherDet);}

  bool CompatibleSig2Noise(const double &Fact) const;
  bool CompatiblePosition(const double &DistMax) const;
  std::string WriteHeader_(ostream & stream, const char*i) const;
  static MatchedDetection* read(istream& r, const char* Format);
  void writen(ostream &s) const ;
};



#include "stringlist.h"

class MatchedDetectionList : public StarList<MatchedDetection>
{
 private:
  StringList imageNames;

 public:
  MatchedDetectionList(const DetectionList &L);
  void ApplyCuts();

  bool OneToOneAssoc(const string &ImageName, DetectionList &L);
};


typedef MatchedDetectionList::iterator MatchedDetectionIterator;
typedef MatchedDetectionList::const_iterator MatchedDetectionCIterator;


#endif /*DETECTION__H */
