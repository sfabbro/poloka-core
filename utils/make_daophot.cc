#include "daophotutils.h"

static void usage(const char *progName)
{
  cerr << " " << progName << " [-m] [-psf] [-cat] <DbImages> \n"
       << "   produces a psf and a fitted catalog with options: \n"
       << "    -psf : produce psf\n"
       << "    -als : produce ALLSTAR catalog\n"
       << "    -man : produce PSF with manual selection of stars\n"
       << "    -c   : combine SExtractor and ALLSTAR catalogs \n"
       << "    -o   : overwrite\n";
}

int main(int nargs, char **args)
{
  if (nargs < 2)  {usage(args[0]);  exit(1);}

  bool overwrite = false;
  bool manual = false;
  bool dopsf = false;
  bool doals = false;
  bool combine = false;

  ReducedImageList toDo;

  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      if ((arg[0] != '-')) 
	{
	  string imName = args[i];
	  ReducedImage *rim = new ReducedImage(imName);
	  if (!rim->ActuallyReduced()) 
	    {
	      cerr << " Image " << imName << " is not fully reduced! \n"; 
	      return EXIT_FAILURE;
	    }
	  toDo.push_back(rim);
	  continue;
	}
      if (strlen(arg) == 2)
	{
	  switch (arg[1])
	    {
	    case 'o' : overwrite = true; break;
	    case 'c' : combine = true; break;
	    default : usage(args[0]); exit(1);
	    }
	  continue;
	}
      if (strlen(arg)> 2)
	{
	  if (strncmp(arg,"-psf",4) == 0) {dopsf = true;continue;}
	  if (strncmp(arg,"-man",4) == 0) {dopsf = true; manual = true; continue;}
	  if (strncmp(arg,"-als",4) == 0) {doals = true;continue;}
	  usage(args[0]); exit(1);
	}      
    }

  if (overwrite && !dopsf && !doals)
    {
      cerr << " DAOPHOT: nothing to overwrite!" << endl;
      usage(args[0]); exit(1);
    }

  for (ReducedImageIterator it=toDo.begin(); it!=toDo.end(); ++it) 
    {
      MakeDaoPsfCat(**it, dopsf, doals, manual, combine, overwrite);
    }

  return EXIT_SUCCESS;
}
