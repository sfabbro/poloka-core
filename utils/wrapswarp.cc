#include <string>
#include <iostream>
#include <fstream>

#include "allreducedimage.h"
#include "swarpstack.h"

static void usage(const char* ProgName)
{
  cerr << ProgName << " [OPTION]... -o DBIMAGE DBIMAGE...\n" 
       << "Wrapper around SWarp\n"
       << "     -r DBIMAGE : specify a astrometric/photometric reference\n"
       << "     -c FILE    : use a SWarp configuration file\n"
       << "     -i FILE    : give a file with a list of DBIMAGEs\n";
  exit(EXIT_FAILURE);
}

int main(int nargs, char **args)
{
  ReducedImageList ril;
  string outName;
  string cardsName;
  ReducedImageRef ref;

  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      if (arg[0] == '-')
	{
	  switch (arg[1]) 
	    {
	    case 'r': 
	      {
		ref = ReducedImageNew(args[++i]); 
		if (ref && !ref->HasImage())
		  {
		    cerr << " the ref image " << ref->Name() << " does not have a fits image\n";
		    return EXIT_FAILURE;
		  }
	      }
	      break;
	    case 'o': outName = args[++i]; break;
	    case 'c': cardsName = args[++i]; break;
	    case 'i': 
	      {
		ifstream ifs(args[++i]);
		if (ifs.is_open()) 
		  {
		    while (ifs.good())
		      {
			string line;
			getline(ifs, line);
			if (line.empty() || line[0] == '#') continue; 
			ReducedImageRef ri = ReducedImageNew(args[i]);
			if (ri->IsValid()) ril.push_back(ri);
			else cerr << " not a valid dbimage: " << arg << endl;
		      }
		    ifs.close();
		  }
	      }
	      break;
	    default: usage(args[0]);
	    }
	  continue;
	}
      else
	{
	  ReducedImageRef ri = ReducedImageNew(args[i]);
	  if (ri->IsValid()) ril.push_back(ri);
	  else cerr << " not a valid dbimage: " << arg << endl;
	}
    }
  

  if (ril.empty())
    {
      cerr << " no input images provided " << endl;
      return EXIT_FAILURE;
    }

  if (outName.empty())
    {
      cerr << " no output name provided " << endl;
      return EXIT_FAILURE;
    }
  
  if (!cardsName.empty()) SetSwarpCardsName(cardsName);
    
  Frame frame;
  if (ref) frame = ref->PhysicalSize();

  SwarpStack ss(outName, ril, ref, frame);
  ss.MakeFits();
  ss.MakeSatur();
  ss.MakeCatalog();

  return EXIT_SUCCESS;
}
