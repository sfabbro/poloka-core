// This may look like C code, but it is really -*- C++ -*-
#ifndef  STARMATCH__H 
#define  STARMATCH__H

#include <iostream>
#include <iterator> // for ostream_iterator
#include <algorithm> // for swap
#include <string>

#include "persistence.h"

#include "point.h"
#include "basestar.h"
#include "gtransfo.h"


/*! \file
   \brief pairs of points

   Used to handles matches of star lists (image/image) or image/catalog
   to fit geometrical and photometric transformations.
   */

//! a pair of stars, usually belonging to different images. 
/*! One would normally assume that
      they are the same object of the sky. The object contains basically two 2d points (called
      later p1 and p2), and
      two pointers (unused by the class and its satellites), that enable in the end to trace back
      the stars in the caller data structures. */

class StarMatch {
  CLASS_VERSION(StarMatch,1);
  #define StarMatch__is__persistent

  friend class StarMatchList;

public: /* if one sets that private, then fitting routines will be a nightmare. we could set all the Transfo classes friend... */
  
  Point point1, point2; //!< 2 points 
  CountedRef<BaseStar> s1, s2; //!< the Star pointers (the pointer is in fact generic, pointed data is never used).
  double distance;

public :
  //! constructor.
  /*! gives 2 points (that contain the geometry), plus pointers to the Star objects
    (which are there for user convenience). */
  StarMatch(const Point &p1, const Point &p2, const BaseStar *S1, const BaseStar *S2) : 
    point1(p1), point2(p2), s1(S1), s2(S2), distance(0.) {};
  
  // the next one would require that StarMatch knows BaseStar which is not mandatory for StarMatch to work
  // StarMatch(BaseStar *S1, BaseStar *S2) : point1(*S1), point2(*S2), s1(S1), s2(S2) {};

  //! returns the distance from T(p1) to p2. 
  double Distance(const Gtransfo &T) const { return  point2.Distance(T.apply(point1)); };

  //! to be used before sorting on distances. 
  void SetDistance(const Gtransfo &T) { distance = Distance(T);};
  //! returns the value computed by the above one.
  double Distance() const { return distance;}

  void Swap() { swap(point1, point2) ; swap(s1,s2) ; }

#ifndef SWIG
  friend bool DecFlux(const StarMatch & S1, const StarMatch & S2)
  {return(S1.s1->flux > S2.s1->flux); }
  friend bool IncreasingDistances(const StarMatch &one, const StarMatch &two) 
  { return(one.distance < two.distance);}
  friend bool DecreasingDistances(const StarMatch &one, const StarMatch &two) 
  { return(one.distance > two.distance);}
#endif
  
  /* comparison that ensures that after a sort, duplicates are next one another */
  explicit StarMatch() {};
  
#ifndef SWIG
  friend ostream& operator << (ostream &stream, const StarMatch &Match); 
#endif

  ~StarMatch() {}


 private:
  //  bool operator <  (const StarMatch & other) const { return (s1 > other.s1) ? s1 > other.s1 : s2 > other.s2;}
 
  /* for unique to remove duplicates */
  bool operator == (const StarMatch &other) const { return (s1 == other.s1 && s2 == other.s2); };
  bool operator != (const StarMatch &other) const { return (s1 != other.s1 || s2 != other.s2); };

  friend bool CompareS1(const StarMatch &one, const StarMatch &two)
    { return ((one.s1 == two.s1) ? (one.distance < two.distance) : ( (const BaseStar*) one.s1 > (const BaseStar *) two.s1));}
  friend bool SameS1(const StarMatch &one, const StarMatch &two)
    { return (one.s1 ==  two.s1);}

  friend bool CompareS2(const StarMatch &one, const StarMatch &two)
    { return ((one.s2 == two.s2) ? (one.distance < two.distance) : ((const BaseStar *)one.s2 > (const BaseStar *) two.s2));}
  friend bool SameS2(const StarMatch &one, const StarMatch &two)
    { return (one.s2 ==  two.s2);}

  //! enables \verbatim cout << mystarMatch << endl; \endverbatim

  //  ClassDef(StarMatch,1);
};

/* =================================== StarMatchList ============================================================ */

#ifndef SWIG
ostream& operator << (ostream &stream, const StarMatch &Match);
#endif


typedef list<StarMatch>::iterator StarMatchIterator;
typedef list<StarMatch>::const_iterator StarMatchCIterator;

//! A list of star matches, 
/*! with the RefineTransfo and Cleanup convenience routines 
         to wrap the fit of the transformation. The fitted transformation
 belongs to the starMatchList object : if you want to store it for later 
use, you have to clone it.
         To choose the type/degree of transfo, see RefineTransfo()
	 */

#ifndef SWIG
ostream& operator << (ostream &stream, const StarMatchList &List);
#endif


class StarMatchList : public list<StarMatch> {


  private :
  int nused;
  int order;
  double chi2;
  Gtransfo *transfo;



  public :


  /* constructor */
  StarMatchList() : nused(0), order(0), chi2(0), transfo(0) {};

  //!carries out a fit with outlier rejection 

  void RefineTransfo(const double &NSigmas, const Gtransfo* PriorTransfo = NULL);

  //! enables to access the fitted transformation. 
  /*! Clone it if you want to store it permanently */
  Gtransfo *Transfo() const { return transfo;}

  //! access to the chi2 of the last fit. 
  double Chi2() const { return chi2;}

  //! computes the chi2 even when there is no fit
  void SetChi2();

  //! returns the number of pairs entering the last fit. 
  int Nused() const { return nused;}
  

  //! returns the degree of freedom for the fit in x and y
  int Dof() const {return (2*nused > transfo->Npar())? (2 *nused-transfo->Npar()): 0;}

  //! returns the order of the used transfo 
  int TransfoOrder() const {return order;}

  
  //! swaps elements 1 and 2 of each starmatch in list. 
    void Swap () ;

  //! returns the average residual. 
  double Residual() const { if (nused != 0) return sqrt(chi2/nused); else return -1;}

  int Cleanup(double DistanceCut, const Gtransfo &Transfo);

  /*! cleans up the list of pairs for pairs that share one of their stars, keeping the closest one.
     The distance is computed using Transfo. */
  int RemoveAmbiguities(const Gtransfo &Transfo);
  
  //! calculates the median of the distance squares of pair lists 
  double MedianDist2() const;
  
  //! sets a transfo between the 2 lists and deletes the previous or default one.  No fit.
  void SetTransfo(const Gtransfo *Transfo) { if (transfo) delete transfo; transfo = Transfo->Clone();}
  //!
  void SetTransfo(const Gtransfo &Transfo) { if (transfo) delete transfo; transfo = Transfo.Clone();}

  //! set transfo according to the given order. 
  void SetTransfoOrder(const int Order);

  /*! returns the inverse transfo (Swap, fit(RefineTransfo) , and Swap). 
         The caller should delete the returned pointer. */
  Gtransfo *InverseTransfo();


  //! Sets the distance (residual) field of all list elements. Mandatory before sorting on distances 
  void SetDistance(const Gtransfo &Transfo);

  //! deletes the tail of the match list 
  void CutTail(const int NKeep);

  //! print the matching transformation quality (transfo, chi2, residual, nused) 
  void DumpTransfo(ostream &stream = cout) const;

  ~StarMatchList() { if (transfo) delete transfo;}
  
 //! write  StarMatchList with a header which  can be read by l2tup 
    void write(const string &filename) const ;
    void write(ostream & pr=cout ) const ;


  //! write  StarMatchList with a header which  can be read by l2tup 
  void write(const Gtransfo &tf, const string &filename) const ;
  void write(const Gtransfo &tf, ostream & pr) const;




    //! enables to dump a match list through : cout << List; 
  void dump(ostream & stream = cout) const { stream << *this; }

  //  ClassDef(StarMatchList,1)

    private:

  StarMatchList(const StarMatchList&); // copies nor properly handled 
  void operator=(const StarMatchList&);


};



#endif /* STARMATCH__H */
