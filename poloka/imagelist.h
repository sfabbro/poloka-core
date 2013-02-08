// This may look like C code, but it is really -*- C++ -*-
#ifndef IMAGELIST__H
#define IMAGELIST__H


#include <list>
#include <iostream>

#include <poloka/countedref.h>

template <class T> class ImageList :
  public  list<CountedRef<T> >
{

private:
  bool shouldDelete;

public:
  //! conctructor. shouldDelete tells if the pointer have to be "deleted" when the list is itself deleted.
  ImageList(const bool ShouldDelete = true) 
      {shouldDelete = ShouldDelete;};

  typedef typename list<CountedRef<T> >::iterator iterator;
  typedef typename list<CountedRef<T> >::const_iterator const_iterator;


  //!
  friend ostream& operator << (ostream & stream, const ImageList& aList) 
    { aList.dump(stream); return stream;}

  //!
  void dump(ostream& stream = cout) const 
  { for (const_iterator ri = this->begin(); ri != this->end(); ++ri) (*ri)->dump(stream);}

};

#include "reducedimage.h"

#endif /* IMAGELIST__H */
