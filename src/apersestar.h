#ifndef APERSESTAR__H
#define APERSESTAR__H

#include <iostream>
#include "sestar.h"

class Image;

#define BAD_GAUSS_MOMENTS 1

// =================================================================
// =================================================================
//
// definition de la classe AperSEStar
//
// =================================================================
// =================================================================
//! A SEstar which has some aperture photometry added
class AperSEStar : public SEStar
{

 public:

  //! gathers data from an aperture photometry measurement
  struct Aperture
  {
    double radius; // in pixels
    double flux; // in ADUS
    double eflux;
    int nbad;
    int ncos; // number of bad pixels flagged as cosmics
    double fcos; // total flux of these pixels

    Aperture() { radius=0; flux=0; eflux=0; nbad=0; ncos=0; fcos=0; }

    bool operator < ( const Aperture &Right) const
    { return radius < Right.radius;}

  };

  vector<Aperture> apers;
  double neighborDist;
  double neighborFlux;
  double gmxx, gmyy, gmxy;
  int gflag;
  
  private :
  void zero() { neighborDist=-1; neighborFlux = -1; gflag=BAD_GAUSS_MOMENTS; gmxx = gmyy = gmxy = 0; }

 public:
  AperSEStar() {zero();}
  AperSEStar(const SEStar &sestar) : SEStar(sestar) {zero();}
  
  //! computes flux (& co) and stores it into an added Aperture instance.
  void ComputeFlux(const Image& I, const Image &W, const Image& C,
		   const double Gain, const double Radius);

  void ComputeShapeAndPos(const Image&I, const Image &W);

  Aperture InterpolateFlux(const double &Radius) const;

  void SetNeighborScores(const BaseStar &Neighbor);
  double NeighborDist() const { return neighborDist;}
  double NeighborFlux() const { return neighborFlux;}


  std::string WriteHeader_(ostream & pr = cout, const char* i = NULL) const;
  void read_it(istream& r, const char* Format);
  static AperSEStar* read(istream& r, const char* Format); 
  void writen(ostream& s) const;
};


ostream& operator << (ostream &stream, const AperSEStar::Aperture &);


#include  "starlist.h"

class AperSEStarList : public StarList<AperSEStar>
{
  
  public:

  AperSEStarList() {};
  
  AperSEStarList(const string &FileName) {read(FileName);}

  explicit AperSEStarList(const SEStarList &L);

    //! routine that sets in each AperSEStar what concerns its closest neighbor.
    /*! important notice: the routine assumes that the objects of
      *this are in the provided List. */
  void SetNeighborScores(const BaseStarList &List, const double Maxdist);

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


bool FindStarShapes(const AperSEStarList &List, double &SizeX, 
		    double &SizeY, double &Corr);


#endif /* APERSESTAR__H */
