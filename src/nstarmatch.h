#ifndef  NSTARMATCH__H 
#define  NSTARMATCH__H
// This may look like C code, but it is really -*- C++ -*-

#include "basestar.h"
#include "gtransfo.h"

/*! May be used to transform coordinates of lists (e.g. for matches)
   without overwriting original ones */
class TStar : public BaseStar
{
  BaseStarRef original;

public:
  TStar(const BaseStar &B, const Gtransfo &T) : BaseStar(B), original(&B) 
  {
    (Point &) *this = T.apply(B);
  }

  TStar(const BaseStar &ActualStar, const BaseStar &JustForPos) : BaseStar(JustForPos), original(&ActualStar) {}
  virtual void writen(ostream &s = cout)const ;
  virtual string WriteHeader_(ostream & stream = cout, const char*i = NULL) const ;
  
  BaseStar const*  get_original() const { return original; }
  BaseStar*  get_original() { return original; }
};
    
class TStarList : public StarList<TStar> {
  public :
    TStarList(const BaseStarList& slist, const Gtransfo& tranfo);
};

typedef StarList<TStar>::iterator TStarIterator;
typedef StarList<TStar>::const_iterator TStarCIterator;

#include <vector>

//! this is hanger of matched objects
class NStarMatch : public BaseStar {

  vector<CountedRef<BaseStar> > stars;
  int nm;

public:
  NStarMatch() : nm(0) {};
  
  //! true size of the match vector
  unsigned int  size() const { return stars.size(); }
  
  //! add one entry into the match at a given index. Position gets updated
  void AddMatch(const BaseStar &Match, const unsigned Index);

  //! removes one entry into the match at a given index. Position gets updated
  void DropMatch(const unsigned Index);


  //! get actual number of matches ( <= size)
  int NumberOfMatches() const { return nm;}


  //! returns the match entered before (if any)
  const BaseStar *GetMatch(const unsigned Index) const;
  BaseStar       *GetMatch(const unsigned Index);


  friend class NStarMatchList;

};

/* ============ NStarMatchList ========================= */


class NStarMatchList : public StarList<NStarMatch> {

 private:

  vector<BaseStarRef> emptyStars;
  vector<string> tags; // output tags for tuple entry names

public :

  int MatchAnotherList(const BaseStarList &List,
		       const double MaxDist, 
		       const BaseStar *EmptyStar,
		       const string &Tag = "");

  int  write(const string &FileName) const;

  int Nlists() const { return int(emptyStars.size());}

  void ApplyCountCut(const int MinCount);

  //! return the tag of the i-th list
  string const&   getTag(int i) const { return tags[i]; }
  
  //! return the tag vector
  vector<string> const& getTags() const { return tags; }

#ifdef STORAGE
  /* constructor */


  //! access to the number of lists. 


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

#endif 
};


typedef NStarMatchList::iterator NStarMatchIterator;
typedef NStarMatchList::const_iterator NStarMatchCIterator;



#endif /* NSTARMATCH__H */
