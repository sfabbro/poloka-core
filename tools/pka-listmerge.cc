#include <list>
#include <string>
#include <iostream>

#include <poloka/dicstar.h>
#include <poloka/reducedimage.h>
#include <poloka/gtransfo.h>
#include <poloka/wcsutils.h>
#include <poloka/fileutils.h>
#include <poloka/fitsimage.h>


static void usage(const char *progname)
{
  cerr << "Usage: " << progname << "[OPTION]... LIST... [-o outfilename] [-t (x->xccd, ra->x)] "
       << "Merge lists\n\n"
       << "    -o FILE : specify output file name (default: out.list)\n"
       << "    -t      : transform coordinates (x->xccd, ra->x...)\n";
  exit(EXIT_FAILURE);
}

int main(int nargs, char **args)
{
  if (nargs<2) usage(args[0]);

  list<string> imList;
  bool transformCoords = false;;
  string outputName="out.list";

  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      if (arg[0] != '-')
	{
	  imList.push_back(arg);
	  continue;
	}
      switch (arg[1])
	{
	case 'o' : i++; outputName = args[i]; break;
	case 't' : transformCoords = true; break;
	default : 
	  cerr << args[0] << ": don't understand " << arg << endl;
	  return EXIT_FAILURE;
	}
    }

  DicStarList outputList;
  bool first =true;
  unsigned firstNKeys = 0;

  bool ok = true;
  for (list<string>::const_iterator it = imList.begin(); it != imList.end(); ++it)
    {
      string listName = *it;
      cout << " considering " << listName << endl;
      string name = DirName(listName);
      ReducedImage ri(name);
      if (!ri.IsValid())
	{
	  cerr << args[0] << ": cannot find reducedimage for " << name << endl;
	  ok = false;
	  continue;
	}
      FitsHeader head(ri.FitsName());
      if (!head.IsValid())
	{
	  cerr << args[0] << ": cannot find fitsimage for  " << name << endl;
	  ok = false;
	  continue;
	}
      int shoot = head.KeyVal("EXPNUM"); // specific to Megacam
      int chip = head.KeyVal("TOADCHIP");
      GtransfoRef wcs = WCSFromHeader(head);
      if (transformCoords && !wcs)
	{
	  cerr << args[0] << ": do not find the expected WCS in " << head.FileName() << endl;
	  ok = false;
	  continue;
	}
      DicStarList starList(listName);
      for (DicStarIterator i = starList.begin(); i != starList.end(); ++i)
	{
	  DicStar &s = **i;
	  if (transformCoords)
	    {
	      s.AddKey("xccd",s.x);
	      s.AddKey("yccd",s.y);
	      // substitute x,y -> ra,dec
	      wcs->apply(s.x,s.y, s.x,s.y);
	    }
	  s.AddKey("chip",chip);
	  s.AddKey("shoot", shoot);

	  if (first)
	    {
	      firstNKeys = s.NKeys();
	      first = false;
	    }
	  else
	    {
	      if (s.NKeys() != firstNKeys)
		{
		  cerr << args[0] << ": no way to merge files with different number of keys\n";
		  return EXIT_FAILURE;
		}
	    }
	    
	  // append to output list 
	  outputList.push_back(*i); 
	}
    }
  cout << " writing " << outputName << endl;
  outputList.write(outputName);
  return ok? EXIT_SUCCESS : EXIT_FAILURE;
}
	      
