#include <iostream>
#include <vector>

#include "psfmatch.h"
#include "fitsimage.h"

// for io
#include "kernelfit_dict.h"
#include "objio.h"
#include "typemgr.h"


static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " best worst  " << endl;
  cerr << "  match PSF of <DbImage(s)> with a reference image (assumed of best seeing) \n" 
       << "   Note:both images have to be on the same pixel grid ! \n";

  exit(-1);
}

int main(int nargs, char **args)
{
  // if nothing is given
  if (nargs < 3){usage(args[0]);}
  
  CountedRef<ReducedImage> refimage = new ReducedImage(args[1]);
  CountedRef<ReducedImage> newimage = new ReducedImage(args[2]);
  PsfMatch psfmatch(*refimage,*newimage);
  psfmatch.FitKernel(true);
  //ReducedImage sub("sub");
  //psfmatch.Subtraction(sub);
  string outputfilename = newimage->Dir()+"/kernel_from_"+refimage->Name()+".xml";
  cout << "image_sub : writing kernel in " << outputfilename << " ..." << endl;
  obj_output<xmlstream> oo(outputfilename);
  oo << *(psfmatch.GetKernelFit());
  oo.close();
  cout << "the end" << endl;
  
  
  return 0;
}
