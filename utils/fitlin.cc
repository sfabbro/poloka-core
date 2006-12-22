
#include <iostream>
#include <string>


#include "imagepsf.h"

static void usage(const string &prog)
{
  cerr << prog << " <dbimage ...> " << endl;
  exit(EXIT_FAILURE);
}
 

int main( int nargs, char **args)
{
  if (nargs <=1) usage(args[0]);
  bool success = true;
  for (int i=1; i < nargs; ++i)
    {
      ReducedImage *ri = new ReducedImage(args[i]);
      if (!ri->IsValid())
	{
	  cout << " cannot find " << args[i] << endl;
	  continue;
	}
      ImagePSF imagePSF(*ri);
      success &=imagePSF.FitNonLinearity();
      //      delete ri; avoid it, already done....
    }

  if (success) return EXIT_SUCCESS; else return EXIT_FAILURE;
}
