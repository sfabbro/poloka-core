#include "psfmatch.h"
#include "kernelfit.h"

static void usage(const char *progName)
{
  cerr << progName << " -r <Ref> <DbImage1> ... <DbImageN> \n"; 
}

int main(int nargs, char **args)
{
  if (nargs <= 4){usage(args[0]); exit(1);}

  ReducedImageList imlist;
  string refName = "NOREF";
  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      if (arg[0] != '-') 
	{
	    string name = string(args[i]);       
	    ReducedImage *im = new ReducedImage(name);
	    imlist.push_back(im);
	    continue;
	}
      switch (arg[1])
	{
	case 'r' :++i; refName = args[i]; break;
	default :{usage(args[0]); exit(1);}
	}
    }	

  if (refName=="NOREF") {usage(args[0]); exit(1);}
  ReducedImageRef ref = new ReducedImage(refName);
  
  for (ReducedImageCIterator it = imlist.begin(); it != imlist.end(); ++it)
    {
      ReducedImage *current = *it;
      PsfMatch match(ref, current);
      match.FitKernel();
      Kernel kern;
      match.KernelToWorst(kern,current->XSize()/2, current->YSize()/2);
      string ker_name = ref->Name()+current->Name()+"_kernel.fits";
      kern.writeFits(ker_name);
    }

  return EXIT_SUCCESS;
}

