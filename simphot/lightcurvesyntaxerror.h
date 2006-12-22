#ifndef LIGHTCURVESYNTAXERROR__H
#define LIGHTCURVESYNTAXERROR__H

#include "polokaexception.h"


class LightCurveSyntaxError: public PolokaException
{
 public:
  LightCurveSyntaxError(const string &Line,const string &Error) : PolokaException("Syntax error")
  {
    append(Error);
    append("in line :\n");
    append(Line);
  }

};

#endif /* LIGHTCURVESYNTAXERROR__H */
