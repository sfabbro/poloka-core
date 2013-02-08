#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include <poloka/dbimage.h>
#include <poloka/reducedimage.h>
#include <poloka/fileutils.h>

static void usage(const char* progname)
{
  cerr << "Usage: " << progname << " [OPTION]... DBIMAGE...\n"
       << "Produce a weight image for DBIMAGE\n\n"
       << "   -o: overwrite\n";
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{ 
  
  ReducedImageList imList;
  bool overwrite = false;
  bool ok = true;

  for (int i=1; i<argc; ++i)
    {
      char *arg = argv[i];
      if ((arg[0] != '-')) 
	{
	  ReducedImage *current = new ReducedImage(argv[i]);
	  if (FileExists(current->FitsName())) 
	    imList.push_back(current);
	  else {
	    cerr << argv[0] << ": " << argv[i] << " not found\n";
	    ok = false;
	  }
	  continue;
	}
      switch (arg[1])
	{
	case 'o' : overwrite = true; break;
	default : usage(argv[0]);
	}
    }
  
  for (ReducedImageIterator it = imList.begin(); it != imList.end(); ++it)
    {
      ReducedImage *current = *it;
      if (overwrite && current->HasWeight())
	ok = remove((current->FitsWeightName()).c_str());
      ok = current->MakeWeight();
    }
  
  return ok? EXIT_SUCCESS : EXIT_FAILURE;
}
