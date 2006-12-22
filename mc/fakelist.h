// This may look like C code, but it is really -*- C++ -*-





#ifndef FAKELIST__H
#define FAKELIST__H

#include <string>
#include <list>
#include <iostream>

#include "countedref.h"
#include "point.h"


using namespace std;

class Frame;

//! lists of Fakes. 
/*!     It is a template class, which means that the Fake
type remains undefined until a user defines it. 
The list related operations (insertion,
sort, traversal) are to be carried out using STL 
list operations. Most of the
Fake operations rely on routines to be provided in 
the Fake class, usually
user defined. The instanciation of this class for 
BaseFake (i.e. the replacement 
of the formal parameter 'Fake' by 'BaseFake') is 
called BaseFakeList. 
Take care: what is stored is pointers on Fake's and 
NOT Fake's. This implies
that Fakes being inserted in the list have to be 
obtained using 'new'. The corresponding
'delete' are invoked in the destructor. */

template<class Fake> class FakeList : public list <CountedRef<Fake> >  {
 
public:
typedef CountedRef<Fake> Element;
typedef typename list<Element>::const_iterator FakeCIterator;
typedef typename list<Element>::iterator FakeIterator;


/* constructors */
//! : default constructor (empty list). 
  FakeList() : list<Element>() {};

 //! reads a FakeList from a file, 
  /*!
     using the read method from the Fake class. 
     See BaseFake for an example of implementation. */
  FakeList(const string &FileName); 
  //! writes to a file
  /*!  calls iteratively the write method of the Fake 
	class. It is unusable if the Fake class does not 
	provide this functionnality.
	see BaseFake to see a possible implementation. 
  not const because the write routines of Root are not*/
  int  write(const string &FileName);

  //! obvious meaning
  int read(const string &FileName);


  void push_back(Fake* t) {list<Element>::push_back(Element(t));} 
  /* the previous one hides the following one ?! */
  void push_back(const Element& e) {list<Element>::push_back(e);}


/* destructor */
  virtual ~FakeList() {};

  
  //! invokes dump(stream) for all Fakes in the list. 
 void dump(ostream &stream = cout ) const { 
    for (FakeCIterator p = this->begin(); 
	 p !=this->end(); ++p) (*p)->dump(stream);}

 
  //! copy the head of the list at the  end of an other list (that may be empty on input)
  void ExtractHead(FakeList<Fake> &Out, int NHead) const;

  //! cuts the end of the list 
  void CutTail(const int NKeep);

  //! copy the part of the list which is included in the frame at the end of another list
 void ExtractInFrame(FakeList<Fake> &Out, const Frame &aFrame) const;
 //! cut the part of the list which is at a distance < mindist of the edges defined by frame.
  void CutEdges(const Frame &aFrame, float mindist);

  //! clears Copy and makes a copy of the list to Copy 
  void CopyTo(FakeList<Fake> &Copy) const;

  //! Clears the list 
  void ClearList() { CutTail(0);};

  //! enables to apply a geometrical transfo if Fake is Basestar or derives from it.
  /*! could be extended to other type of transformations. */

  template<class Operator> void ApplyTransfo(const Operator &Op) 
    {for (FakeIterator p = this->begin(); p != this->end(); ++p) (*p)->Apply(Op);}

  //! returns the closest Fake from a given location. 
  Fake* FindClosest(double X, double Y) const; 
  Fake* FindClosest(double X, double Y,FakeList<Fake> const  & stlparallel, 
		    Fake* & star2  ) const ;

  //! returns the closest Fake from a given location. 
  //! and removes it from the given list.
  Fake* FindAndRemoveClosest(double X, double Y);

  //! same as above. Can be used with any of our star-like stuff.
  Fake* FindClosest(const Point &P) const { return FindClosest(P.x, P.y);};

#ifndef __CINT__
  //! true if location has a nearby star in a ring between mindist and maxdist
 bool HasCloseNeighbor(double X, double Y, double maxdist, double mindist=0.5) const;

  //! same as above. Can be used with any of our star-like stuff.
  bool HasCloseNeighbor(const Point &P, double maxdist, double mindist=0.1) const
   {return HasCloseNeighbor(P.x, P.y, maxdist, mindist);}

  //! nearby star to a star but not itself
  Fake* ClosestNeighbor(double X, double Y, double mindist=0.1) const;

  //! same as above. Can be used with any of our star-like stuff.
  Fake* ClosestNeighbor(const Point &P, double mindist=0.1) const
   {return ClosestNeighbor(P.x, P.y, mindist);}

  int NumberOfNeighbors(const double &X, const double &Y, const double &distmax) const;
  int AllNeighbors(FakeList &NeighborList, const double &X, const double &Y, const double &distmax) const;
  int AllNeighbors(FakeList &NeighborList, const Point &Pt, const double &distmax) const
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

  private :
  // no implicit copies.
#ifndef SWIG
  // FakeList(const FakeList<Fake>&);
  // FakeList& operator = (const FakeList&);
  // Authorized IMPLICIT COPIES !! so that when cloning imagesubtraction, for MonteCarlo simulation, it can copies the KernelFit, which contains a FakeList.(DH)
#endif

};

  //! enables \verbatim  cout << my_list; \endverbatim
template <class Fake>  ostream & operator <<(ostream &stream, const FakeList<Fake> &List) 
    {List.dump(stream); return stream; };



#ifdef USE_ROOT
#include "rootstuff.h"
template <class Fake> class FakeListWithRoot : public FakeList<Fake>, public TObject
{
public :

  FakeListWithRoot();
  FakeListWithRoot(const string &FileName);
  FakeListWithRoot(const char* FileName) {this->read(string(FileName));}
  int write( const string &FileName);
  int read(const string &FileName);
  
private :
  int root_read(const string &FileName);
  int root_write(const string &FileName);

  //  bool root_write(const string &FileName);
  ClassDefT(FakeListWithRoot,1)
    };

#endif /* USE_ROOT */  

#endif /* FAKELIST__H */
