#ifndef DBCONFIGEXCEPTION__H
#define DBCONFIGXCEPTION__H
 
#include <iostream>
#include <poloka/polokaexception.h>


/* in polokaexception.h, you may read about the simple rules
   we assume for exception handling.
*/
 
//! class specialized in DbConfig exceptions (syntax errors in the dbconfig file) 
class DbConfigException : public PolokaException
{
 public:

  DbConfigException(const string &Mess) : PolokaException("(DBCONFIG) : "+Mess) {}

};

#endif /* DBCONFIGEXCEPTION__H */
