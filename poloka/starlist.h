// This may look like C code, but it is really -*- C++ -*-
#ifndef STARLIST__H
#define STARLIST__H

#include <string>
#include <list>
#include <iostream>

#include <poloka/countedref.h>
#include <poloka/point.h>
#include <poloka/globalval.h>


using namespace std;

class Frame;
class fastifstream;

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
  GlobalVal glob;
 
public:
typedef CountedRef<Star> Element;
typedef typename list<Element>::const_iterator StarCIterator;
typedef typename list<Element>::iterator StarIterator;


/* constructors */
//! : default constructor (empty list). 
  StarList() {};

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
  */

  Star *EmptyStar() const { return new Star();}

  int  write(const string &FileName) const;

  //! obvious meaning
  int read(const string &FileName);

  //! enables to access global values (lines starting with '@' in ascii files)
  GlobalVal &GlobVal() { return glob;}

  //! enables to access global values (lines starting with '@' in ascii files)
  const GlobalVal &GlobVal() const { return glob;}

  /* the previous one hides the following one ?! */
  void push_back(const Element& e) {list<Element>::push_back(e);}


/* destructor */
  virtual ~StarList() {};

  
  //! invokes dump(stream) for all Stars in the list. 
 void dump(ostream &stream = cout ) const { 
    for (StarCIterator p = this->begin(); 
	 p !=this->end(); ++p) (*p)->dump(stream);}

  //!a model routine to sort the list 
  /*! see DecreasingFlux() to see what it is, if you 
     want another sorting criterion) */
  // le premier de la liste a le plus grand flux
  void FluxSort();

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
  {for (StarIterator p = this->begin(); p != this->end(); ++p) Op.TransformStar(*(*p));}

  //! returns the closest Star from a given location. 
  Star* FindClosest(double X, double Y) const;

  //! same as above. Can be used with any of our star-like stuff.
  Star* FindClosest(const Point &P) const { return FindClosest(P.x, P.y);};

  //! returns the closest Star from a given location. 
  //! and removes it from the given list.
  Star* FindAndRemoveClosest(double X, double Y);

#ifndef __CINT__
  //! true if location has a nearby star in a ring between mindist and maxdist
  //! if minflux>0, then the condition flux > minflux is required.
  bool HasCloseNeighbor(double X, double Y, double maxdist, double mindist=0.1,double minflux=-1) const;

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

  int read(fastifstream & rd);


public :

  //!: with descriptor for l2tup 
  int  write(ostream & pr) const;


};

  //! enables \verbatim  cout << my_list; \endverbatim
#ifndef SWIG
template <class Star>  ostream & operator <<(ostream &stream, const StarList<Star> &List) 
    {List.dump(stream); return stream; }
#endif 



#endif /* STARLIST__H */
