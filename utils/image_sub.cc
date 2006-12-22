#include <iostream>
#include <fstream>
#include <vector>

#include "psfmatch.h"
#include "fitsimage.h"
#include "imagesubtraction.h"

// for io
//#include "kernelfit_dict.h"
//#include "objio.h"
//#include "typemgr.h"


static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " best worst  " << endl;
  cerr << "options (to put after best worst) : -s  do the subtraction" << endl;
  cerr << "  match PSF of <DbImage(s)> with a reference image (assumed of best seeing) \n" 
       << "   Note:both images have to be on the same pixel grid ! \n";

  exit(-1);
}

int main(int nargs, char **args)
{
  // if nothing is given
  if (nargs < 3){usage(args[0]);}
  
  bool makesub = false;
  
  CountedRef<ReducedImage> refimage = new ReducedImage(args[1]);
  CountedRef<ReducedImage> newimage = new ReducedImage(args[2]);

  for (int i=3; i<nargs; ++i) {
    char *arg = args[i];
    if (arg[0] != '-')
      usage(args[0]);
    switch (arg[1])
      {
      case 's' : makesub = true; break;
      default : usage(args[0]);
      }
  }


  // write kernel in xml file
  string outputfilename = newimage->Dir()+"/kernel_from_"+refimage->Name()+".dat";
  cout << "image_sub : writing kernel in " << outputfilename << " ..." << endl;
  //  obj_output<xmlstream> oo(outputfilename);
  ofstream oo(outputfilename.c_str());
  double phoratio = 1;
  double chi2 = 0.;
  int nstars = 0;
  int nparams = 0;
  
  if(makesub) {
    cout << "image_sub : creating sub " << refimage->Name() << " " << newimage->Name() << endl;
    ImageSubtraction sub("sub",refimage,newimage);
    sub.MakeFits();
    cout << "image_sub : writing kernel in " << outputfilename << " ..." << endl; 
    sub.GetKernelFit()->write(oo);
    //    oo <<*(sub.GetKernelFit());

    phoratio = sub.PhotomRatio();
    chi2 = sub.Chi2();
    nstars = sub.Nstars();
    nparams = sub.Nparams();
    
  }else{
    cout << "image_sub : psfmatch "  << refimage->Name() << " " << newimage->Name() << endl;
    PsfMatch psfmatch(refimage,newimage,NULL,true);
    if( ! psfmatch.FitKernel(true) ) {
      cout << "psfmatch failed" << endl;
      
      // dump something for log
      return EXIT_FAILURE;
    }
    cout << "image_sub : writing kernel in " << outputfilename << " ..." << endl; 
    //    obj_output<xmlstream> oo(outputfilename);
    psfmatch.GetKernelFit()->write(oo);
    //    oo << *(psfmatch.GetKernelFit());

    phoratio = psfmatch.PhotomRatio();
    chi2 = psfmatch.Chi2();
    nstars = psfmatch.Nstars();
    nparams = psfmatch.Nparams();
  }
  oo.close();
  
  FitsHeader head(newimage->FitsName(), RW);
  head.AddOrModKey("PMRATIO", phoratio);
  head.AddOrModKey("PMCHI2" , chi2);
  head.AddOrModKey("PMNSTAR", nstars);
  head.AddOrModKey("PMNPAR" , nparams);

  cout << "the end" << endl;
  
  
  return EXIT_SUCCESS;
}
