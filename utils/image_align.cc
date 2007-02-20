#include <iostream>
#include <vector>

#include "transformedimage.h"
#include "fitsimage.h"
#include "polokaexception.h"

static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " [OPTIONS] <DbImage(s)> " << endl;
  cerr << "  register and resample images relatively to a geometric reference image\n" 
       << "    OPTIONS: \n"
       << "    -geo <name>: indicate the geometric reference image <name>. Default is first one \n"
       << "    -no_catalog: no transformed catalog \n"
       << "    -no_weight: no transformed weight\n"
       << "    -no_satur: no transformed satur\n"
       << "    -wcs     : use wcs for alignement (instead of image catalogs)\n" ;
   
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
  
  int todo = DoFits | DoCatalog | DoSatur | DoWeight;
  bool use_wcs = false;
  
  
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
      if (strcmp(arg,"no_catalog")==0) {todo&=(!DoCatalog); continue;}
      if (strcmp(arg,"no_weight")==0) {todo&=(!DoWeight); continue;}
      if (strcmp(arg,"no_satur")==0) {todo&=(!DoSatur); continue;}
      if (strcmp(arg,"wcs")==0) {use_wcs=true; continue;}
      
      // unrecognized option
      usage(args[0]);      
    }
  
  if (geoName == "NOGEO") geoName = toAlign.front()->Name();
  try {
    ReducedImage geoRef(geoName);
    ReducedImageList aligned;
  
    unsigned int n = ImagesAlign(toAlign, geoRef, aligned, todo, use_wcs);
    if(n!=toAlign.size()) {
      throw(PolokaException("Not all images were aligned"));
    }

  }catch(PolokaException p)
    {
      p.PrintMessage(cout);
      return EXIT_FAILURE;	      
    } 
  
  return EXIT_SUCCESS;
}

