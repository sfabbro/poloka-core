// This may look like C code, but it is really -*- C++ -*-



#ifndef STARLIST__H
#define STARLIST__H

#include <string>
#include <list>
#include <iostream>

#include "countedref.h"
#include "point.h"     // NRL: 04/2004


using namespace std;

class Frame;

//! lists of Stars. 
/*!     It is a template class, which means that the Star
type remains undefined until a user defines it. 
The list related operations (insertion,
sort, traversal) are to be carried out using STL 
list operations. Most of the
Star operations rely on routines to be provided in 
the Star class, usually
user defined. The instanciation of this class for 
BaseStar (i.e. the replacement 
of the formal parameter 'Star' by 'BaseStar') is 
called BaseStarList. 
Take care: what is stored is pointers on Star's and 
NOT Star's. This implies
that Stars being inserted in the list have to be 
obtained using 'new'. The corresponding
'delete' are invoked in the destructor. */

template<class Star> class StarList : public list <CountedRef<Star> >  {
 
public:
typedef CountedRef<Star> Element;
typedef typename list<Element>::const_iterator StarCIterator;
typedef typename list<Element>::iterator StarIterator;


/* constructors */
//! : default constructor (empty list). 
  StarList() : list<Element>() {};

 //! reads a StarList from a file, 
  /*!
     using the read method from the Star class. 
     See BaseStar for an example of implementation. */
  StarList(const string &FileName); 
  //! writes to a file
  /*!  calls iteratively the write method of the Star 
	class. It is unusable if the Star class does not 
	provide this functionnality.
	see BaseStar to see a possible implementation. 
  not const because the write routines of Root are not*/
  int  write(const string &FileName);

  //! obvious meaning
  int read(const string &FileName);


  void push_back(Star* t) {list<Element>::push_back(Element(t));} 
  /* the previous one hides the following one ?! */
  void push_back(const Element& e) {list<Element>::push_back(e);}


/* destructor */
  virtual ~StarList() {};

  
  //! invokes dump(stream) for all Stars in the list. 
 void dump(ostream &stream = cout ) const { 
    for (StarCIterator p = begin(); 
	 p !=end(); ++p) (*p)->dump(stream);}

  //!a model routine to sort the list 
  /*! see DecreasingFlux() to see what it is, if you 
     want another sorting criterion) */
  // le premier de la liste a le plus grand flux
  void FluxSort() { sort(&DecreasingFlux);}

  //! copy the head of the list at the  end of an other list (that may be empty on input)
  void ExtractHead(StarList<Star> &Out, int NHead) const;

  //! cuts the end of the list 
  void CutTail(const int NKeep);

  //! copy the part of the list which is included in the frame at the end of another list
  void ExtractInFrame(StarList<Star> &Out, const Frame &aFrame) const;
 //! cut the part of the list which is at a distance < mindist of the edges defined by frame.
  void CutEdges(const Frame &aFrame, float mindist);

  //! clears Copy and makes a copy of the list to Copy 
  void CopyTo(StarList<Star> &Copy) const;

  //! Clears the list 
  void ClearList() { CutTail(0);};

  //! enables to apply a geometrical transfo if Star is Basestar or derives from it.
  /*! could be extended to other type of transformations. */

  template<class Operator> void ApplyTransfo(const Operator &Op) 
    {for (StarIterator p = begin(); p !=end(); ++p) (*p)->Apply(Op);}

  //! returns the closest Star from a given location. 
  Star* FindClosest(double X, double Y) const;

  //! same as above. Can be used with any of our star-like stuff.
  Star* FindClosest(const Point &P) const { return FindClosest(P.x, P.y);};

#ifndef __CINT__
  //! true if location has a nearby star in a ring between mindist and maxdist
  bool HasCloseNeighbor(double X, double Y, double maxdist, double mindist=0.1) const;

  //! same as above. Can be used with any of our star-like stuff.
  bool HasCloseNeighbor(const Point &P, double maxdist, double mindist=0.1) const
   {return HasCloseNeighbor(P.x, P.y, maxdist, mindist);}

  //! nearby star to a star but not itself
  Star* ClosestNeighbor(double X, double Y, double mindist=0.1) const;

  //! same as above. Can be used with any of our star-like stuff.
  Star* ClosestNeighbor(const Point &P, double mindist=0.1) const
   {return ClosestNeighbor(P.x, P.y, mindist);}

  int NumberOfNeighbors(const double &X, const double &Y, const double &distmax) const;
  int AllNeighbors(StarList &NeighborList, const double &X, const double &Y, const double &distmax) const;
  int AllNeighbors(StarList &NeighborList, const Point &Pt, const double &distmax) const
  {return AllNeighbors(NeighborList, Pt.x, Pt.y, distmax);}
  int NumberOfNeighbors(const Point &Pt, const double &distmax) const
  {return NumberOfNeighbors(Pt.x, Pt.y, distmax);}
#endif


protected :
  int ascii_read(const string &FileName);

private :

  int read(istream & rd);

  //!: with descriptor for l2tup 
  int  write(ostream & pr) const;

  //!: without descriptor for l2tup 
  int  write_wnoheader(ostream & pr) const;

public :

  // NRL 02/2004 We need it for SWIG...
  //  private :
  // no implicit copies.
#ifndef SWIG
     StarList(const StarList<Star>&);
     StarList& operator = (const StarList&);
#endif

};

  //! enables \verbatim  cout << my_list; \endverbatim
template <class Star>  ostream & operator <<(ostream &stream, const StarList<Star> &List) 
    {List.dump(stream); return stream; };



#ifdef USE_ROOT
#include "rootstuff.h"
template <class Star> class StarListWithRoot : public StarList<Star>, public TObject
{
public :

  StarListWithRoot();
  StarListWithRoot(const string &FileName);
  StarListWithRoot(const char* FileName) {this->read(string(FileName));}
  int write( const string &FileName);
  int read(const string &FileName);
  
private :
  int root_read(const string &FileName);
  int root_write(const string &FileName);

  //  bool root_write(const string &FileName);
  ClassDefT(StarListWithRoot,1)
    };

#endif /* USE_ROOT */  

#endif /* STARLIST__H */
