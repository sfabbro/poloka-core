#ifndef JIMSTAR_SEEN
#define JIMSTAR_SEEN

#include <iostream>
#include <fstream>

#include "basestar.h"
//#include "astroutils.h"
//#include "sestar.h"

/*! \file */
//! Jim Calibration Star 
/*! To match an image with the  catalogue of Jim 
 * format:  ra dec i r g z cgal  (cgal=0 ==> etoile)
*/

class JimStar : public BaseStar {
  public :
  
  
    JimStar();
  JimStar(double xx, double yy, double ff); 
  
  string Name() const {return name;};
  double Ra() const {return ra;};
  double Dec() const {return dec;};
  
  double Imag() const {return imag;};
  double Rmag() const {return rmag;};
  double Gmag() const {return gmag;};
  double Zmag() const {return zmag;};
  double Cgal() const {return cgal;};

  //! to read once the object is created 
  virtual void    Read(istream& r, const char *Format); 
  
  //! to read and create the object  
  static JimStar* read(istream& r, const char *Format); 
  
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
  double ra; // R.A. in deg
  double dec; // Dec in deg
  double imag;
  double rmag;
  double gmag;
  double zmag;
  double cgal;

#ifdef USE_ROOT
  ClassDef(JimStar,1) // Root's stuff
#endif 
    };
    
    
/********************   FIN DEFINITION JimStar   **********************/


/* what concerns the JimStarList's : */
#include "starlist.h"
#include "fitsimage.h"

//@{
//! definition of the list and iterators
#ifdef USE_ROOT
typedef StarListWithRoot<JimStar> JimStarList;
#else
typedef StarList<JimStar> JimStarList;
#endif /* USE_ROOT */

typedef StarList<JimStar>::const_iterator JimStarCIterator;
typedef StarList<JimStar>::iterator JimStarIterator;


JimStarList *GetSelectedJimStarList(const FitsHeader &header, const Frame &W,const string &standardfile);

#endif
