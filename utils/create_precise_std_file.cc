#include <math.h>

#include "standardstar.h"
#include "usnoutils.h"
#include "wcsutils.h"

#ifndef M_PI
#define  M_PI           3.14159265358979323846  /* pi */
#endif


int main(int argc, char **argv)
{

  char *env_var = getenv("STANDARDFILE");
  if (!env_var)
    {
      cerr << " you should define STANDARDFILE env var to run this code " << endl;
      return -1;
    }
  string standardfile(env_var);

  StandardStarList stdstarlist(standardfile);

  if (stdstarlist.size()==0)
    {
      cout << "Bad standard file !!!" << endl;
      return -1;
    }

  for (StandardStarIterator si = stdstarlist.begin(); si != stdstarlist.end(); ++si)
    {
      StandardStar * pstar = (StandardStar *) *si;
      double RaDeg = pstar->Ra();
      double DecDeg= pstar->Dec();
      double RaMin = RaDeg-0.5; 
      double RaMax = RaDeg+0.5;
      double DecMin= DecDeg-0.5;
      double DecMax= DecDeg+0.5;
      
      BaseStarList usnoLocalStarList;
      UsnoRead(RaMin,RaMax,DecMin,DecMax,RColor,usnoLocalStarList);
      ConvertMagToFlux(&usnoLocalStarList);
      usnoLocalStarList.FluxSort();

      BaseStar *usnoStar = usnoLocalStarList.front();

      pstar->x = usnoStar->x;
      pstar->y = usnoStar->y;
    }

  string newname = standardfile + "_new";
  stdstarlist.write(newname);
}
