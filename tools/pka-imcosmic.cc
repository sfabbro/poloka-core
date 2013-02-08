#include <poloka/reducedimage.h>
#include <poloka/fitsimage.h>
#include <poloka/fileutils.h>
#include <poloka/sestar.h>
#include <poloka/seeing_box.h>


static void usage(const char* progname)
{
  cerr << "Usage: " << progname << " [OPTION]... DBIMAGE...\n"
       << "Produce a cosmic mask for DBIMAGE\n\n"
       << "   -o: overwrite\n";
  exit(EXIT_FAILURE);
}

int main(int nargs, char ** args)
{
  if (nargs<2) usage(args[0]);

  ReducedImageList imList;
  bool overwrite = false;
  bool ok = true;

  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      if (arg[0] != '-')
	{
	  ReducedImage *current = new ReducedImage(args[i]);
	  if (FileExists(current->FitsName())) 
	    imList.push_back(current);
	  else {
	    cerr << args[0] << ": " << args[i] << " not found\n";
	    ok = false;
	  }
	  continue;
	}
      switch (arg[1])
	{
	case 'o' : overwrite = true; break;
	default : usage(args[0]);
	}
    }
  
  for(ReducedImageIterator it = imList.begin(); it != imList.end(); ++it)
    {
      ReducedImage *current = *it;
      if (overwrite && current->HasCosmic())
	ok = remove((current->FitsCosmicName()).c_str());
      ok = current->MakeCosmic();
    }
  
  return ok? EXIT_SUCCESS : EXIT_FAILURE;
}
