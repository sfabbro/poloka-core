#include <iostream>
#include <cstdio> // for remove
#include <string>

#include <vector>

#include "reducedimage.h"

static void usage(const char *prog)
{
  cerr << prog << " [-o (overwrite)] <dbimages> " << endl
       << " makes the perture photometry catalogue of dbimage(s) " << endl;
  exit (-1);
}
  


int main(int nargs, char **args)
{
  vector<string> names;
  bool overwrite = false;

  for (int i=1; i< nargs; i++)
    {
      char *arg = args[i];
      if (arg[0] != '-') 
	{
	  names.push_back(arg);
	  continue;
	}
      switch (arg[1])
	{
	case 'h' :
	  usage(args[0]);
	  break;
	case 'o' :
	  overwrite = true ;
	  break;
	default:
	  cerr << " don't understand " << arg << endl;
	  usage(args[0]);
	}
      
    } // end loop on arguments

  bool ok = true;
  for (unsigned k=0; k < names.size(); ++k)
    {
      ReducedImage ri(names[k]);
      if (!ri.IsValid())
	{
	  ok = false;
	  cerr << " could not find DbImage " << names[k] << endl;
	  continue;
	}
      if (overwrite)
	{
	  string fileName(ri.AperCatalogName());
	  remove(fileName.c_str());
	}
      ok &= ri.MakeAperCat();
    }
  return (ok=true) ? EXIT_SUCCESS : EXIT_FAILURE;
}
