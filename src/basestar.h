
// This may look like C code, but it is really -*- C++ -*-


#ifndef BASESTAR__H
#define BASESTAR__H


#include <iostream>
#include <cstdio>
#include <string>


#include "point.h"
#include "countedref.h"

using namespace std;


// AVANT :
// Le centre du pixel (i=0,j=0) (coordonnees IJ) 
// a pour coordonnees  (x=0.5, y=0.5) (coordonnees XY) 
// ainsi une etoile en (x=0, y=0) est centree ds le coin en bas a gauche du
//
// APRES : (i=0,j=0) == (x=0., y=0.) pour tout le monde
// !! SAOimage (i=0,j=0) correspond a (1,1)


#define MEMPIX2DISK 1


#include "rootstuff.h" // empty include when not using Root

#define DECALAGE_IJ_XY 0.
#define DECALAGE_XY_IJ 0.
/*! \file */

//! The base class for handling stars. Used by all matching routines.
class BaseStar : public Point, public RefCount
{
  CLASS_VERSION(BaseStar,1);
  #define BaseStar__is__persistent
  public : // si quelqu'un connait un moyen efficace d'eviter ca...
double flux; 


  public:
 BaseStar(){x=0;y=0;flux=0;};
  //! constructor
 BaseStar(double xx, double yy, double ff) : Point(xx,yy), flux(ff) 
  {};
 BaseStar(const Point &a_point, double a_flux) : Point(a_point), flux(a_flux)
  {};

  //! access stuff.
  double X() const { return x;}
  //! .
  double Y() const { return y;}

  static BaseStar* read(istream & rd, const char *format);

  void  read_it(istream & rd, const char *format);
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


  //  ClassDef(BaseStar,1) // no ";" ....
};




//! enables to sort easily a starList (of anything that derives from BaseStar) 
bool DecreasingFlux(const BaseStar *S1, const BaseStar *S2);



int DecodeFormat(const char *FormatLine, const char *StarName);


/* what conscerns the BaseStarList's : */
#include "starlist.h"

typedef StarList<BaseStar> BaseStarList;



typedef BaseStarList::const_iterator BaseStarCIterator;
typedef BaseStarList::iterator BaseStarIterator;
typedef CountedRef<BaseStar> BaseStarRef;

#endif /* BASESTAR__H */  
