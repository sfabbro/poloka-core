// This may look like C code, but it is really -*- C++ -*-
#ifndef IMAGELIST__H
#define IMAGELIST__H


#include <vector>
#include <iostream>

#include "countedref.h"

#include "rootstuff.h"


template <class T> class ImageList : public  vector< CountedRef<T> > 
#ifdef USE_ROOT
, public TObject 
#endif /* USE_ROOT */
{

private:
  bool shouldDelete;

public:
  //! conctructor. shouldDelete tells if the pointer have to be "deleted" when the list is itself deleted.
  ImageList(const bool ShouldDelete = true) 
      {shouldDelete = ShouldDelete;};

  typedef typename vector<CountedRef<T> >::iterator iterator;
  typedef typename vector<CountedRef<T> >::const_iterator const_iterator;


  //!
  friend ostream& operator << (ostream & stream, const ImageList& aList) 
    { aList.dump(stream); return stream;}

  //!
  void dump(ostream& stream = cout) const 
  { for (const_iterator ri = begin(); ri != end(); ++ri) (*ri)->dump(stream);}

#ifndef SWIG
  ClassDef(ImageList,1);
#endif /* SWIG */

};

#include "reducedimage.h"

#endif /* IMAGELIST__H */
