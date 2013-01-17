#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "allreducedimage.h"
#include "swarpstack.h"

static void usage(const char* ProgName)
{
  cerr << ProgName << " [OPTION]... -o DBIMAGE DBIMAGE...\n" 
       << "Wrapper around SWarp\n"
       << "     -s         : stack subtraction instead of calibrated images\n"
       << "     -o DBIMAGE : output stack DBIMAGE\n"
       << "     -r DBIMAGE : specify a astrometric/photometric reference\n"
       << "     -c FILE    : specify a SWarp configuration file instead of default\n"
       << "     -i FILE    : file with a list of DBIMAGEs instead of argument list\n\n";
  exit(EXIT_FAILURE);
}

int main(int nargs, char **args)
{

  if (nargs <= 1) usage(args[0]);

  ReducedImageList ril;
  string outName, cardsName;
  ReducedImageRef ref;
  map<string,Point> dcrval;
  DbImageKind imageType = Calibrated;

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
	    case 's': imageType = Subtracted; break;
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
			string name;
			istringstream iline(line);
			iline >> name;
			ReducedImageRef ri = ReducedImageNew(name);
			if (ri->IsValid()) ril.push_back(ri);
			else cerr << " not a valid dbimage: " << arg << endl;
			Point pt;
			if (iline >> pt.x >> pt.y)
			  dcrval[name] = pt;
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
  ss.SetSwarpType(imageType);
  ss.dcrval = dcrval;
  ss.MakeFits();
  ss.MakeSatur();
  ss.MakeCatalog();

  return EXIT_SUCCESS;
}
