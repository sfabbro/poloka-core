// This may look like C code, but it is really -*- C++ -*-

#include <iostream>
#include <iterator>
#include <list>
#include <string>

#include "point.h"
#include "basestar.h"
#include "gtransfo.h"
#ifndef  NSTARMATCH__H 
#define  NSTARMATCH__H


template <class T> class Table2D
{

private :
  int nx;
  int ny;
  T* data;

public :

  Table2D(const int Nx,const int Ny):nx(Nx), ny(Ny){data= new T[nx*ny];};
  Table2D(const int Nx,const int Ny, T init_val):nx(Nx), ny(Ny)
  {data= new T[nx*ny]; for (int i=0; i<nx*ny; ++i) data[i] = init_val;};
  T operator () (const int i, const int j) const {return data[i+j*nx];};
  T& operator() (const int i, const int j) {return data[i+j*nx];};
  ~Table2D() {delete [] data;};
  
};



class NStarMatch : public BaseStar {

friend class NStarMatchList;

  public :

  // Because we try to match n lists and because not necessary all stars will have a representant 
  // in each list, we need 2 numbers :
  int npointers;
  int actualMatches; 
  // int ref;
  Point* points;
  const BaseStar** star;
//! The distance will be calculated from observation i to the first observation.
  public :

  NStarMatch(int n);

 
  NStarMatch(const BaseStar *Star, const int i, const int n);

  bool StarExist(const int i) const {if (!star[i]) return false; return true;}

  double x(const int i) const {return star[i]->x;}

  double y(const int i) const {return star[i]->y;}

  double flux(const int i) const {return star[i]->flux;}

  double Distance(const int i, const int j, const Gtransfo &T) const 
  { if (star[i] && star[j]) return  points[j].Distance(T.apply(points[i]));
  else return 1e8;};

  // these 2 routines become useless if C++ vectors are used instead of
  // "manual" arrays
  ~NStarMatch();
  NStarMatch(const NStarMatch &O);

  // private:


};

/* =================================== NStarMatchList ============================================================ */

typedef list<NStarMatch>::iterator NStarMatchIterator;
typedef list<NStarMatch>::const_iterator NStarMatchCIterator;

class NStarMatchList : public list<NStarMatch> {


private :
  int nlists;
  Table2D<int> nused;
  Table2D<int> order;
  Table2D<double> chi2;
  Table2D<Gtransfo*> transfo;

  BaseStarList AllStars;
  NStarMatchList(const NStarMatchList&); // to forbid copies
  const BaseStar **emptyStar;

public :


  //  void operator=(const StarMatchList&);

  /* constructor */
  NStarMatchList(const int n);

#ifdef STORAGE
  NStarMatchList(int n, BaseStarList &List) : nlists(n),nused(n,n),order(n,n),chi2(n,n),transfo(n,n,NULL)
  {
    for (BaseStarCIterator si = List.begin(); si != List.end(); ++si)
      {
	NStarMatch nstrm(*si,0,n);
	push_back(nstrm);
      }
  }
#endif 

  //! enables to access the fitted transformation. 
  /*! Clone it if you want to store it permanently */
  Gtransfo* Transfo(int i, int j) const { return transfo(i,j);}
  Gtransfo*& Transfo(int i, int j) { return transfo(i,j);}

  //! access to the number of lists. 
  int Nlists() const { return nlists;}

  //! access to the chi2 of the last fit. 
  double Chi2(int i, int j) const { return chi2(i,j);}

  //! returns the number of pairs entering the last fit. 
  int Nused(int i, int j) const { return nused(i,j);}
  
  //! returns the order of the used transfo 
  int TransfoOrder(int i, int j) const {return order(i,j);}

  //! returns the average residual. 
  double Residual(int i, int j) const { if (nused(i,j) != 0) return sqrt(chi2(i,j)/nused(i,j)); else return -1;}

  //! sets a transfo between the 2 lists and deletes the previous or default one.  No fit.
  void SetTransfo(const Gtransfo *Transf, const int i, const int j) 
  { 
    //    if (transfo(i,j)) delete transfo(i,j); 
    transfo(i,j) = Transf->Clone();
  }
  //!
  void SetTransfo(const Gtransfo &Transf, const int i, const int j) 
  { 
    //    if (transfo(i,j)) delete transfo(i,j); 
    transfo(i,j) = Transf.Clone();
  }

  Point ApplyTransfo(const int i, const int j, Point pt);
  
  /*! EmptyStar is only used by write. It has to be of the same type as
    the stars in List. The default constructor is usually OK. */
  void MatchAnotherList(const BaseStarList &List, const int rank, const double MaxDist, const BaseStar *EmptyStar );
  void MatchAnotherList2(const BaseStarList &List, const int rank, const double MaxDist, const BaseStar *EmptyStar );
  
  ~NStarMatchList()
  { 
    for (int i=0; i<nlists ; i++) for (int j=0; j<nlists ; j++) if (transfo(i,j)) delete transfo(i,j);
    for (int i=0; i<nlists; ++i) if (emptyStar[i]) delete emptyStar[i];
  delete [] emptyStar;
  }
  
 //! write  StarMatchList with a header which  can be read by l2tup 
    void short_write(const string &filename) const ;
    void short_write(ostream & pr=cout ) const ;

 //! write  StarMatchList with a header which  can be read by l2tup 
    void write(const string &filename) const ;
    void write(ostream & pr=cout ) const ;

    //! enables to dump a match list through : cout << List; 
  friend ostream& operator << (ostream &stream, const NStarMatchList &List)
      { stream << " number of elements " << List.size() << endl; copy(List.begin(), List.end(), ostream_iterator<NStarMatch>(stream)); return stream;} 

  void dump(ostream & stream = cout) const { stream << cout; }
};

#endif /* NSTARMATCH__H */
