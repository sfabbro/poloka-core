#ifndef STANDARDSTAR_SEEN
#define STANDARDSTAR_SEEN

#include <iostream>
#include <fstream>

#include "basestar.h"
#include "astroutils.h"
#include "sestar.h"

// To match an image with the Standard catalogue


typedef enum StandardColor {NONE, VBAND, BBAND, UBAND, RBAND, IBAND};


class StandardStar : public BaseStar {

public :

  double airmass; 
  double fluxpersec;
  double efluxpersec;
  int flag;
  StandardColor color;

  StandardStar();
  StandardStar(double xx, double yy, double ff); 

  string Name() const {return name;};
  double AirMass() const {return airmass;};
  double FluxPerSec() const {return fluxpersec;};
  int Flag() const {return flag;};
  double Ra() const {return RaStringToDeg(ra);};
  double Dec() const {return DecStringToDeg(dec);};

  double BVmag() const {return bvmag;};
  double UBmag() const {return ubmag;};
  double VRmag() const {return vrmag;};
  double RImag() const {return rimag;};
  double VImag() const {return vimag;};


  double Vmag() const {return vmag;};
  double Bmag() const {return vmag + bvmag;};
  double Umag() const {return vmag + bvmag + ubmag;};
  double Rmag() const {return vmag - vrmag;};
  double Imag() const {return vmag - vimag;};

  double dVmag() const {return dvmag;};
  double dBmag() const {return sqrt(dvmag*dvmag + dbvmag*dbvmag);};
  double dUmag() const {return sqrt(dvmag*dvmag + dbvmag*dbvmag + dubmag*dubmag);};
  double dRmag() const {return sqrt(dvmag*dvmag + dvrmag*dvrmag);};
  double dImag() const {return sqrt(dvmag*dvmag + dvimag*dvimag);};

  double Magnitude(const StandardColor couleur);
  double Magnitude();
  double eMagnitude(const StandardColor couleur);
  double eMagnitude();
  double DeltaColor();

   //! to read once the object is created 
  virtual void    Read(istream& r, const char *Format); 

   //! to read and create the object  

  static StandardStar* read(istream& r, const char *Format); 

   //!  to write the StarList header with the string appended to every ntuple variable (with no end)
  string WriteHeader_(ostream & pr = cout, const char *i=NULL) const;
 
   //! for dump with NO end-of-line
  virtual void    dumpn(ostream& s = cout) const;

   //! for dump  
  virtual void    dump(ostream& s = cout) const ;
  
   //! for write with NO end-of-line
  virtual void    writen(ostream& s = cout) const ;

   //! for write 
  virtual void    write(ostream& s= cout)  const ;

private :

void Set_to_Zero();

protected :
  string name;
  string ra;
  string dec;
  double vmag;
  double bvmag;
  double ubmag;
  double vrmag;
  double rimag;
  double vimag;
  int n,m;
  double dvmag;
  double dbvmag;
  double dubmag;
  double dvrmag;
  double drimag;
  double dvimag;

#ifndef SWIG
  ClassDef(StandardStar,1) // Root's stuff
#endif

};


/********************   FIN DEFINITION StandardStar   **********************/


/* what concerns the StandardStarList's : */
#include "starlist.h"
#include "fitsimage.h"

//@{
//! definition of the list and iterators
#ifdef USE_ROOT
typedef StarListWithRoot<StandardStar> StandardStarList;
#else
typedef StarList<StandardStar> StandardStarList;
#endif /* USE_ROOT */

typedef StarList<StandardStar>::const_iterator StandardStarCIterator;
typedef StarList<StandardStar>::iterator StandardStarIterator;

#ifndef SWIG
BaseStarList* Standard2Base(StandardStarList * This);
const BaseStarList* Standard2Base(const StandardStarList * This);
#endif


StandardColor GetColor (const FitsHeader &header);
StandardStarList *GetSelectedStandardStarList(const FitsHeader &header);
int GetStandardZeroPoint(StandardStarList *standardList, SEStarList &sestarlist, const FitsHeader &header, double *zeropoint, double *err_zp, int &nzero);
double FitZeroPoint(StandardStarList &stdstarlist,double &a1, double &a2, double mean, double sig, double &err_zp, double &err_a1, double &err_a2, string outfilename);
double ReadTheoZeroPoint(const FitsHeader &header); 

#endif
