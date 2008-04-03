#ifndef MULTIMAGSTAR__H
#define MULTIMAGSTAR__H

#include <iostream>


#include "sestar.h"
#include "apersestar.h"
#include "ellipticaper.h"


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
    double seeing ;
    double f_aper;
    double ef_aper ;
    double f_aper_other ;
    ShortMagBox(){f_circ = 0. ; ef_circ = 0. ;f_auto = 0. ; ef_auto = 0. ; m_auto = 0. ; em_auto = 0. ;  seeing=0. ;f_aper =0.;  ef_aper=0. ;  f_aper_other=0. ;}
};


class MultiMagSEStar : public SEStar {


 public:
    // pour conserver info d'association a un autre catalogue (par ex : cfhtls)
    CountedRef<BaseStar> star ;
    double star_dist ;


  vector<ShortMagBox> magboxes;
  double alpha ;
  double delta ;
  double x_orig ;
  double y_orig ;

  // ouverture de photometrie sur l'image i
  Elliptic_Aperture ell_aper ;
  Elliptic_Aperture g_ell_aper ;

  // parametre de forme sur l'image de detction.
  double gx ;
  double gy ;
  double gmxx ;
  double gmyy ;
  double gmxy ;
  double gmxx_loc ;
  double gmyy_loc ;
  double gmxy_loc ;

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
bool UpDate_Assoc(SEStarList &L, string band) ;
 bool UpDate(SEStarList &L, string band);
 void ComputeMag(int kbox, string band, double ZP, double eZP);
 MultiMagSEStar*  FindEllipticClosest(double xx, double yy, double dilatation, double RadMin, double Radius) const ;

 int GetBandNumber(string band) ;
 void SetBandNumber(string band, int n);
 int GetNBand() ;
 void SetNBand(int n);

  //! the standard one plus some checks
 void check() const ;
 int write(const std::string &FileName) const;

 int  MatchToOtherList(BaseStarList *l);
 
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
