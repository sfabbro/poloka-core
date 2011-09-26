#include "objecttofit.h"

#include "lightcurvesyntaxerror.h"
#include "fileutils.h" // DecomposeString et al

ObjectToFit::ObjectToFit(const string& Line)
{  
  jdmin = -1.e+30 ; jdmax= 1.e+30 ;
  vector<string> words;
  DecomposeString(words, Line);
  if (words.size() < 2)
    throw(LightCurveSyntaxError(Line,"not enough words"));


  type=0; // default is star+galaxy
  x = atof(words[0].c_str());
  y = atof(words[1].c_str());
  vector<string>::const_iterator it=words.begin();
  ++it; ++it; // skip the 2 coordinates
  for ( ; it != words.end(); ++it)
    {
      vector<string> option;
      DecomposeString(option, *it, "=");
      if      (option[0]=="DATE_MIN") jdmin = atof(option[1].c_str());
      else if (option[0]=="DATE_MAX") jdmax = atof(option[1].c_str());
      else if (option[0]=="NAME")     name  = option[1];
      else if (option[0]=="TYPE")     type  = atoi(option[1].c_str());
      else if (option[0]=="BAND")     band  = option[1][0];
      
      else 
	throw(LightCurveSyntaxError(Line, " unknown argument : "+option[0]));
    }

  if ((jdmin > 0) && (jdmax > 0) && (jdmin >= jdmax))
    throw(LightCurveSyntaxError(Line, " bad dates "));
}

