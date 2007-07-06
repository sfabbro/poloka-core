#ifndef MULTIMAGSTAR__H
#define MULTIMAGSTAR__H

#include <iostream>


#include "sestar.h"
#include "apersestar.h"

using namespace std;

class Image;
class FitsHeader;

class fastifstream ;


class CalibBox
{
public:
    string band ;
  double ZP;
  double sigZP;
  double ZPell;
  double sigZPell;
  double seeing;
 
  CalibBox() { band=""; ZP=0. ; sigZP =0. ; ZPell=0. ; sigZPell =0. ; seeing =0. ;}
  
  void SetParameters(const string Band, const double zeropoint, const double sigmaZP, const double Seeing) { band=Band; ZP = zeropoint; sigZP = sigmaZP ; seeing = Seeing ;}

  };



class ShortMagBox
{
public:
    CalibBox calib ;
    double f_auto ;
    double ef_auto ;
    double f_circ ;
    double ef_circ ;
    double m_auto ;
    double em_auto ;
    ShortMagBox(){f_circ = 0. ; ef_circ = 0. ;f_auto = 0. ; ef_auto = 0. ; m_auto = 0. ; em_auto = 0. ;}
};


class MultiMagSEStar : public SEStar {


 public:
  vector<ShortMagBox> magboxes;
  double alpha ;
  double delta ;
  double x_orig ;
  double y_orig ;

  // variable utilisee pour les sort des liste de galaxies voisines d'une SN
  //pas dans les routines write et read pour l'instant
  double ell_dist ;
  double norm_dist ;
  double dist ;

 public:
  MultiMagSEStar(){SetToZero();}
  MultiMagSEStar(const SEStar &sestar);
 
 
 void ComputeMag(int kbox, string band, double ZP, double eZP);
 double SqEllipticDistance(double xx, double yy, double dilatation, double RadMin,double Radius) const;
 double NormalizedDistance(double xx, double yy, double dilatation, double RadMin,double Radius) const;

  std::string WriteHeader_(ostream & pr = cout, const char* i = NULL) const;
  void read_it(fastifstream& r, const char* Format);
  static MultiMagSEStar* read(fastifstream& r, const char* Format); 
  void writen(ostream& s) const;
  //! for dump with NO end-of-line
  virtual void    dumpn(ostream& s = cout) const;

 private :
  void SetToZero();
};


bool 
IncSqEllipticDist(const MultiMagSEStar *S1, const MultiMagSEStar *S2);
bool 
IncNormalizedDistance(const MultiMagSEStar *S1, const MultiMagSEStar *S2);

#include  "starlist.h"

class MultiMagSEStarList : public StarList<MultiMagSEStar>
{
  
  public:

  MultiMagSEStarList() {}
  
  MultiMagSEStarList(const string &FileName); 

  explicit MultiMagSEStarList(const SEStarList &L);

 void CopySEStarList(const SEStarList &L);
 void ComputeAlphaDelta(const FitsHeader & head);

 bool UpDate(SEStarList &L, string *bands, int nband);
 void ComputeMag(int kbox, string band, double ZP, double eZP);
 MultiMagSEStar*  FindEllipticClosest(double xx, double yy, double dilatation, double RadMin, double Radius) const ;



  //! the standard one plus some checks
 void check() const ;
 int write(const std::string &FileName) const;
 
};


// pour determiner l'hote du SN : distance elliptique 
// rappel : parametrisation de l'ellipse de l'objet :  cxx x^2 + cxy xy + cyy y^2 = 1
// correspond a une ellipse de grad axe A() pixels etc.
//  d^2 = cxx x^2 + cxy xy + cyy y^2. 
//!! si on a choisi un cercle, cxx=cyy=1, cxy=0 et d est alors la distance en pixels !

void FindEllipticNeighb(double xsn, double ysn, double dist,
			MultiMagSEStarList & stlin, 
			MultiMagSEStarList & stl_neighb);
// pour determiner l'hote du SN : distance elliptique normalisee
// l'ellipse dephotometrie est definie par  cxx x^2 + cxy xy + cyy y^2 = R^2
// avec R = kron factor = 2.5 * kron_radius, si < 3.5, alors = 3.5
// et si (kron_factor * sqrt(A*B) <= Radmin * 0.5 (Radmin=8 pix ici)  alors on prend 1 cercle de rayon 
// ici fixe a Radius=8 pix.
// d_norm = d/R si ellipse, =d/Radius si cercle.

void FindNormalizedDistNeighb(double xsn, double ysn, double dist,
			MultiMagSEStarList & stlin, 
			MultiMagSEStarList & stl_neighb);
// routine de check : quand les 2 distances donnent des resultats differents.
// a eliminer plus tard.
void CheckNeighbFinders(string name, double xsn, double ysn, double dist,
			MultiMagSEStarList & stlin);


const BaseStarList* MultiMagSE2Base(const MultiMagSEStarList * This);

BaseStarList& MultiMagSE2Base(MultiMagSEStarList &This);

const BaseStarList& MultiMagSE2Base(const MultiMagSEStarList &This);

SEStarList* MultiMagSE2SE(MultiMagSEStarList * This);

const SEStarList* MultiMagSE2SE(const MultiMagSEStarList * This);

SEStarList& MultiMagSE2SE(MultiMagSEStarList &This);

const SEStarList& MultiMagSE2SE(const MultiMagSEStarList &This);

/* ================================== */

AperSEStarList* MultiMagSE2AperSE(MultiMagSEStarList * This);
const AperSEStarList* MultiMagSE2AperSE(const MultiMagSEStarList * This);
AperSEStarList& MultiMagSE2AperSE(MultiMagSEStarList &This);
const AperSEStarList& MultiMagSE2AperSE(const MultiMagSEStarList &This);




typedef MultiMagSEStarList::const_iterator MultiMagSEStarCIterator;
typedef MultiMagSEStarList::iterator MultiMagSEStarIterator;
typedef CountedRef<MultiMagSEStar> MultiMagSEStarRef;



#endif /* MULTIMAGSTAR__H */
