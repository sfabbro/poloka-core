#include "daophot.h"
#include "daophotutils.h"

static void usage(const char *progName)
{
  cerr << " " << progName << " [OPTIONS] <DbImages> \n"
       << "   produces a psf and a fitted catalog with options: \n"
       << "    -psf : produce psf\n"
       << "    -als : produce ALLSTAR catalog\n"
       << "    -m   : merge ALLSTAR catalog into se.list \n"
       << "    -o   : overwrite\n";
}

int main(int nargs, char **args)
{
  if (nargs < 2)  {usage(args[0]);  exit(1);}

  bool over  = false;
  bool dopsf = false;
  bool doals = false;
  bool merge = false;

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
	    case 'o' : over  = true; break;
	    case 'm' : merge = true; break;
	    default : usage(args[0]); exit(1);
	    }
	  continue;
	}
      if (strlen(arg)> 2)
	{
	  if (strncmp(arg,"-psf",4) == 0) { dopsf = true; continue; }
	  if (strncmp(arg,"-als",4) == 0) { doals = true; continue; }
	  usage(args[0]); exit(1);
	}      
    }

  if (over && !dopsf && !doals)
    {
      cerr << " " << args[0] << ": nothing to overwrite!" << endl;
      usage(args[0]); exit(1);
    }


  if ((dopsf) && (!doals))
    for (ReducedImageIterator it=toDo.begin(); it!=toDo.end(); ++it) 
      //MakeDaoPsf(**it, over);
      MakeExperimentalPsf(**it);
  else if (doals) 
    for (ReducedImageIterator it=toDo.begin(); it!=toDo.end(); ++it) 
      MakeDaoPsfAls(**it, merge, over);

  return EXIT_SUCCESS;
}
