#include <iostream>
#include <vector>

#include "transformedimage.h"
#include "fitsimage.h"

static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " [OPTIONS] <DbImage(s)> " << endl;
  cerr << "  register and resample images relatively to a geometric reference image\n" 
       << "    OPTIONS: \n"
       << "    -geo <name>: indicate the geometric reference image <name>. Default is first one \n";
  exit(-1);
}

int main(int nargs, char **args)
{
  // check defaults
  // if nothing is given
  if (nargs < 2){usage(args[0]);}
  if (nargs == 2){cerr << " Align at least 2 images !!\n\n"; usage(args[0]);}
  ReducedImageList toAlign;
  string geoName("NOGEO");
  
  // loop over arguments
  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];

      // images
      if (arg[0] != '-') 
	{
	  ReducedImage *im = new ReducedImage(arg);
	  toAlign.push_back(im);
	  continue;
	}

      // options
      arg++;
      if (strcmp(arg,"geo")==0) { ++i; geoName = args[i]; continue;}

      // unrecognized option
      usage(args[0]);      
    }

  if (geoName == "NOGEO") geoName = toAlign.front()->Name();
  ReducedImage geoRef(geoName);
  ReducedImageList aligned;
  ImagesAlign(toAlign, geoRef, aligned, DoFits | DoCatalog | DoSatur | DoWeight);

  return EXIT_SUCCESS;
}

