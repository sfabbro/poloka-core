#include "dicstar.h"
#include "reducedimage.h"
#include "gtransfo.h"
#include "wcsutils.h"
#include "fileutils.h" // DirName
#include "fitsimage.h"

#include <vector>
#include <string>
#include <iostream>

using namespace std;

static void usage(const char *prog)
{
  cerr << prog << " <lists> [-o outfilename] [-t (x->xccd, ra->x)] " << endl;
  exit(-1);
}

int main(int nargs, char **args)
{
  if (nargs<2) usage(args[0]);
  vector<string> listNames;
  bool transformCoords = false;;
  string outputName="out.list";
  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      if (arg[0] != '-')
	{
	  listNames.push_back(arg);
	  continue;
	}
      switch (arg[1])
	{
	case 'o' : i++; outputName = args[i]; break;
	case 't' : transformCoords = true; break;
	default : cerr << " don't understand " << arg << endl; 
	  usage(args[0]); break;
	}
    }

  DicStarList outputList;
  bool first =true;
  unsigned firstNKeys = 0;
  for (unsigned k=0; k < listNames.size(); ++k)
    {
      string &listName = listNames[k];
      cout << " considering " << listName << endl;
      string name = DirName(listName);
      ReducedImage ri(name);
      if (!ri.IsValid())
	{
	  cout << " cannot find reducedimage for  " << name << endl;
	  continue;
	}
      FitsHeader head(ri.FitsName());
      if (!head.IsValid())
	{
	  cout << " cannot find fitsimage for  " << name << endl;
	  continue;
	}
      int shoot = head.KeyVal("EXPNUM"); // specific to Megacam
      int chip = head.KeyVal("TOADCHIP");
      Gtransfo *wcs = NULL;
      if (transformCoords && !WCSFromHeader(head, wcs))
	{
	  cerr << " do not find the expected WCS in " << head.FileName() << endl;
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
		  cout << " no way to merge files with different number of keys !!! " << endl;
		  cout << " stopping here " << endl;
		  return EXIT_FAILURE;
		}
	    }
	    
	  // append to output list 
	  outputList.push_back(*i); 
	}
    }
  cout << " writing " << outputName << endl;
  outputList.write(outputName);
  return EXIT_SUCCESS;
}
	      
