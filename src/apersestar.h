#ifndef APERSESTAR__H
#define APERSESTAR__H

#include <iostream>
#include "sestar.h"

class Image;

#define BAD_GAUSS_MOMENTS 1
#define BAD_GAUSS_POS_VARIANCE 2




class fastifstream;

//
// definition of Aperture
//

//! gathers data from an aperture photometry measurement
struct Aperture
{
    double radius; // in pixels
    double flux; // in ADUS
    double eflux;
    int nbad;
    int ncos; // number of bad pixels flagged as cosmics
    double fcos; // total flux of these pixels
    double fother; // flux of other objects in aperture
    Aperture() { radius=0; flux=0; eflux=0; nbad=0; ncos=0; fcos=0; fother=0;}
    void computeflux(double x, double y, 
		     const Image& I, const Image& W, const Image *pC, const Image *pS, 
		     const double Gain, const double Radius,int segmentation_number);
    bool operator < ( const Aperture &Right) const
    { return radius < Right.radius;}

};


//
// definition of AperSEStar
//

//! A SEstar which has some aperture photometry added
class AperSEStar : public SEStar
{

 public:
  vector<Aperture> apers;
  double neighborDist;
  double neighborFlux;
  bool   neighborContamination;
  double maxFluxContamination;
  double gmxx, gmyy, gmxy;
  int gflag;
  
  private :
  void zero() { 
    neighborDist=-1; 
    neighborFlux = -1; 
    neighborContamination=false; 
    maxFluxContamination=-1; 
    gflag=BAD_GAUSS_MOMENTS; gmxx = gmyy = gmxy = 0; }


 public:  


void computeflux(const Image& I, const Image& W, const Image *C, 
		 const Image *S, 
		 const double Gain, const double Radius, bool sort_radii);


 AperSEStar() {zero();}
 AperSEStar(const SEStar &sestar) : SEStar(sestar) {zero();}
 // ~AperSEStar() { __NB_APERSESTARS__ -= 1; } 
 
  //! computes flux (& co) and stores it into an added Aperture instance.
      // sort_radii = true: apertures are sorted in increasing radii order
  // 2 images : image + weight image
  // 4 images : image + weight image + cosmic image + segmentation image
  void ComputeFlux(const Image& I, const Image &W, const double Gain, const double Radius, bool sort_radii=true);
  // to match the old routine
  void ComputeFlux(const Image& I, const Image &W, const Image& Cosmic, 
		   const Image &Seg, const double Gain, const double Radius,
		   bool sort_radii=true);
  void ComputePos(const Image&I, const Image &W);
  void ComputeShapeAndPos(const Image&I, const Image &W, 
			  const double &Gain);

  Aperture InterpolateFlux(const double &Radius) const;

  void SetNeighborScores(const BaseStar &Neighbor);
  double NeighborDist() const { return neighborDist;}
  double NeighborFlux() const { return neighborFlux;}


  std::string WriteHeader_(ostream & pr = cout, const char* i = NULL) const;
  void read_it(fastifstream& r, const char* Format);
  static AperSEStar* read(fastifstream& r, const char* Format); 
  void writen(ostream& s) const;
};


ostream& operator << (ostream &stream, const Aperture &);


#include  "starlist.h"

class AperSEStarList : public StarList<AperSEStar>
{
  
  public:

  AperSEStarList() {};
  
  explicit AperSEStarList(const string &FileName) {read(FileName);}

  explicit AperSEStarList(const SEStarList &L);

    //! routine that sets in each AperSEStar what concerns its closest neighbor.
    /*! important notice: the routine assumes that the objects of
      *this are in the provided List. */
  void SetNeighborScores(const BaseStarList &List, const double Maxdist);
  int NSuccess();

  //! the standard one plus some checks
  int write(const std::string &FileName) const;


 
};


BaseStarList* AperSE2Base(AperSEStarList * This);

const BaseStarList* AperSE2Base(const AperSEStarList * This);

BaseStarList& AperSE2Base(AperSEStarList &This);

const BaseStarList& AperSE2Base(const AperSEStarList &This);

SEStarList* AperSE2SE(AperSEStarList * This);

const SEStarList* AperSE2SE(const AperSEStarList * This);

SEStarList& AperSE2SE(AperSEStarList &This);

const SEStarList& AperSE2SE(const AperSEStarList &This);




typedef AperSEStarList::const_iterator AperSEStarCIterator;
typedef AperSEStarList::iterator AperSEStarIterator;
typedef CountedRef<AperSEStar> AperSEStarRef;


class StarScoresList;

double  Find_Fluxmax_Min(const AperSEStarList &List,  const double MinSN, const double frac_elim, const double histo_val_min, int verbose);


bool FindStarShapes(const AperSEStarList &List, double MinSN, double &SizeX, 
		    double &SizeY, double &Corr, StarScoresList &Scores, 
		    const double frac_elim=0., const double histo_val_min=0., int verbose=0);


#endif /* APERSESTAR__H */
