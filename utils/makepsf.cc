
#include <iostream>
#include <string>


#include "imagepsf.h"
#include "toadscards.h"
#include "polokaexception.h"

static void usage(const string &prog)
{
  cerr << prog << " [-f (force)] [-L (use local star cat)] [-d <datacards>] <dbimage ...> " << endl;
  exit(EXIT_FAILURE);
}
 

int main( int nargs, char **args)
{
  if (nargs <=1) usage(args[0]);
  bool success = true;
  vector<string> names;
  bool force = false, use_external_cat= true;
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
	case 'L' : use_external_cat= false; break;
	case 'd' : SetDatacardsFileName(args[++i]); break;
	default : usage(args[0]); break;
	}
    }
  try
    {
  for (unsigned k=0; k<names.size(); ++k)
    {
      success &= MakePSF(names[k], force, use_external_cat);
    }
    }
  catch(PolokaException p)
    {
      p.PrintMessage(cout);
      success = false;
    };

 
  if (success) return EXIT_SUCCESS; else return EXIT_FAILURE;
}
