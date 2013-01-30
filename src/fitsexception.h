#ifndef FITSEXCEPTION__H
#define FITSEXCEPTION__H
 
#include <iostream>
#include "polokaexception.h"


/* in polokaexception.h, you may read about the simple rules
   we assume for exception handling.
*/
 
//! class specialized in FITS exceptions (namely cfitsio errors) 
class FitsException : public PolokaException
{
 public:

  FitsException(const string &Mess) : PolokaException("FITS: "+Mess) {}

};

#endif /* FITSEXCEPTION__H */
