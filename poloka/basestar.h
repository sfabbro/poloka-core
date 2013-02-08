
// This may look like C code, but it is really -*- C++ -*-


#ifndef BASESTAR__H
#define BASESTAR__H


#include <iostream>
#include <cstdio>
#include <string>

#include <poloka/fatpoint.h>
#include <poloka/countedref.h>
#include <poloka/starlist.h>


using namespace std;


class fastifstream;

// AVANT :
// Le centre du pixel (i=0,j=0) (coordonnees IJ) 
// a pour coordonnees  (x=0.5, y=0.5) (coordonnees XY) 
// ainsi une etoile en (x=0, y=0) est centree ds le coin en bas a gauche du
//
// APRES : (i=0,j=0) == (x=0., y=0.) pour tout le monde
// !! SAOimage (i=0,j=0) correspond a (1,1)


#define MEMPIX2DISK 1


#define DECALAGE_IJ_XY 0.
#define DECALAGE_XY_IJ 0.
/*! \file */

// tell other pieces of code that BaseStar now derives from FatPoint (instead of Point)
#define BASESTAR_HAS_POSITION_ERRORS

//! The base class for handling stars. Used by all matching routines.
class BaseStar : public FatPoint, public RefCount
{

  public : // si quelqu'un connait un moyen efficace d'eviter ca...
  double flux,eflux; 


  public:
  BaseStar(){x=0;y=0;flux=0; eflux=0;}
  //! constructor
  BaseStar(const double& xx, const double& yy, const double& ff, const double& eff=0.) 
  : FatPoint(xx,yy), flux(ff), eflux(eff) {}
  BaseStar(const Point &a_point, const double&ff, const double& eff=0.) 
    : FatPoint(a_point), flux(ff), eflux(eff) {}

  //! access stuff.
  double X() const { return x;}
  //! .
  double Y() const { return y;}

  static BaseStar* read(fastifstream & rd, const char *format);

  void  read_it(fastifstream & rd, const char *format);
  virtual void write(ostream &s = cout)const ;
  virtual void writen(ostream &s = cout)const ;

#ifndef SWIG
  //! allows cout << aBaseStar;
  friend ostream& operator << (ostream &stream, const BaseStar &s)
  { s.dump(stream); return stream;}
#endif

  virtual void dump(ostream & stream = cout) const { stream << "x "<< x << " y " << y << " flux " << flux << endl;}
  virtual void dumpn(ostream & stream = cout) const { stream << "x "<< x << " y " << y << " flux " << flux << " ";}

#ifndef SWIG
  BaseStar& operator=(const Point &P) {this->x = P.x; this->y = P.y; return (*this);};
#endif

  static const char *TypeName() { return "BaseStar";}

  virtual ~BaseStar(){};

  virtual string WriteHeader_(ostream & stream = cout, const char*i = NULL) const ;

  virtual void WriteHeader(ostream & stream = cout) const;
#ifdef USE_ROOT

#endif
  //  ClassDef(BaseStar,1) // no ";" ....
};

//! Number of values read for this format
unsigned NValsBaseStar(const char *Format);


//! enables to sort easily a starList (of anything that derives from BaseStar) 
bool DecreasingFlux(const BaseStar *S1, const BaseStar *S2);



int DecodeFormat(const char *FormatLine, const char *StarName);


/* what conscerns the BaseStarList's : */

typedef StarList<BaseStar> BaseStarList;

typedef BaseStarList::const_iterator BaseStarCIterator;
typedef BaseStarList::iterator BaseStarIterator;
typedef CountedRef<BaseStar> BaseStarRef;

#endif /* BASESTAR__H */  
