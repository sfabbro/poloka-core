#ifndef MATCHEXCEPTION__H
#define MATCHEXCEPTION__H
 
#include <iostream>
#include "polokaexception.h"


/* in polokaexception.h, you may read about the simple rules
   we assume for exception handling.
*/
 
//! class specialized in MATCH exceptions (namely cmatchio errors) 
class MatchException : public PolokaException
{
 public:

  MatchException(const string &Mess) : PolokaException("MATCH: "+Mess) {}

};

#endif /* MATCHEXCEPTION__H */
