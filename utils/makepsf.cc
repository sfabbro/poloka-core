
#include <iostream>
#include <string>


#include "imagepsf.h"

static void usage(const string &prog)
{
  cerr << prog << " [-f (force)] <dbimage ...> " << endl;
  exit(EXIT_FAILURE);
}
 

int main( int nargs, char **args)
{
  if (nargs <=1) usage(args[0]);
  bool success = true;
  vector<string> names;
  bool force = false;
  for (int i=1; i < nargs; ++i)
    {
      const char *arg = args[i];
      if (arg[0] != '-')
	{
	  names.push_back(args[i]);
	  continue;
	}
      switch (arg[1])
	{
	case 'h' : usage(args[0]); break;
	case 'f' : force = true; break;
	default : usage(args[0]); break;
	}
    }
  for (unsigned k=0; k<names.size(); ++k)
    {
      success &= MakePSF(names[k], force);
    }
  if (success) return EXIT_SUCCESS; else return EXIT_FAILURE;
}
