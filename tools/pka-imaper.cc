#include <iostream>
#include <cstdio>
#include <string>
#include <list>

#include <poloka/reducedimage.h>
#include <poloka/polokaexception.h>

static void usage(const char *progname)
{
  cerr << "Usage: " << progname << "[-o] DBIMAGE...\n"
       << "Performs aperture photometry for DBIMAGE\n\n"
       << "    -o : overwrite\n\n";
  exit(EXIT_FAILURE);
}

int main(int nargs, char **args)
{
  if (nargs<2) usage(args[0]);
  list<string> imList;
  bool overwrite = false;

  for (int i=1; i< nargs; i++)
    {
      char *arg = args[i];
      if (arg[0] != '-') 
	{
	  imList.push_back(arg);
	  continue;
	}
      switch (arg[1])
	{
	case 'o' :
	  overwrite = true;
	  break;
	default:
	  cerr << args[0] << ": don't understand " << arg << endl;
	  return EXIT_FAILURE;
	}
      
    } // end loop on arguments

  bool ok = true;
  for (list<string>::const_iterator it = imList.begin(); it != imList.end(); ++it)
    {
      try {
	string name = *it;
	ReducedImage ri(name);
	if (!ri.IsValid())
	  {
	    ok = false;
	    cerr << args[0] << ": could not find DbImage " << name << endl;
	    continue;
	  }
	if (overwrite)
	  {
	    ok = remove(ri.AperCatalogName().c_str());
	    ok = remove(ri.StarCatalogName().c_str());
	  }
      ok &= ( ri.MakeAperCat() && ri.MakeStarCat() );
      
      } catch(PolokaException p) {
	p.PrintMessage(cout);
	ok = false;
      }
    }
  return ok? EXIT_SUCCESS : EXIT_FAILURE;
}
