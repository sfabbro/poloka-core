#include <string>
#include <iostream>

#include "swarpstack.h"

static void usage(const std::string &ProgName)
{
  std::cerr << ProgName << " -o <outName> [-r <RefName>] [-c cardsName] dbimages..." 
	    << std::endl;
  std::cerr << "          when a ref is provided, the output image is geometrically aligned with ref " << std::endl;
  exit(1);
}

int main(int nargs, char **args)
{
  ReducedImageList ril;
  std::string outName;
  std::string cardsName;
  ReducedImage *ref = NULL;
  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      if (arg[0] == '-')
	{
	  switch (arg[1]) 
	    {
	    case 'r' : ++i; ref = new ReducedImage(args[i]); break;
	    case 'o' : ++i; outName = args[i]; break;
	    case 'c' : ++i; cardsName = args[i]; break;
	    default : usage(args[0]);
	    }
	  continue;
	}
      else
	{
	  ReducedImage *ri = new ReducedImage(args[i]);
	  if (ri->IsValid()) ril.push_back(ri);
	  else
	    {
	      std::cerr <<" cannot find " << args[i] << std::endl;
	      delete ri;
	    }
	}
    }

  if (ref && !ref->HasImage())
    {
      std::cerr << " the ref image " << ref->Name() << " does not have a fits image !!" 
		<< std::endl;
      exit(1);
    }
  if (ril.size() == 0)
    {
      std::cerr << " no input images provided " << std::endl;
      usage(args[0]);
    }
  if (outName == "")
    {
      std::cerr << " No output name provided " << std::endl;
      usage(args[0]);
    }

  if (cardsName  != "") SetSwarpCardsName(cardsName);
    

  Frame frame;
  if (ref) frame = ref->PhysicalSize();

  SwarpStack ss(outName, ril, ref, frame);
  ss.MakeFits();
  ss.MakeSatur();
  ss.MakeCatalog();
  return EXIT_SUCCESS;
}
