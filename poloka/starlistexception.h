#ifndef STARLISTEXCEPTION__H
#define STARLISTEXCEPTION__H
 
#include <iostream>
#include <poloka/polokaexception.h>


/* in polokaexception.h, you may read about the simple rules
   we assume for exception handling.
*/
 
//! class specialized in MATCH exceptions (namely cmatchio errors) 
class StarListException : public PolokaException
{
 public:

  StarListException(const string &Mess) : PolokaException("StarList: "+Mess) {}

};

#endif /* STARLISTEXCEPTION__H */
