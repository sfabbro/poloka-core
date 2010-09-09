#include <iostream>
#include <fstream>
#include <vector>

#include "fitsimage.h"
#include "imagesubtraction.h"
#include "polokaexception.h"
#include "toadscards.h"
#include "kernelfit.h"

// for io
//#include "kernelfit_dict.h"
//#include "objio.h"
//#include "typemgr.h"


static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " best worst  " << endl;
  cerr << "options (to put after best worst) : [-s]  do the subtraction" << endl;
  cerr << "[-c <datacardsName>] : use these datacards " << endl;

  cerr << "  match PSF of <DbImage(s)> with a reference image (assumed of best seeing) \n" 
       << "   Note:both images have to be on the same pixel grid ! \n";

  exit(-1);
}

int main(int nargs, char **args)
{
  // if nothing is given
  if (nargs < 3){usage(args[0]);}
  string dataCardsName;
  
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
      case 'c' : ++i; dataCardsName = args[i]; break;
      default : usage(args[0]);
      }
  }

  try {

    if (dataCardsName != "") SetDatacardsFileName(dataCardsName);
    // write kernel in xml file
    string outputfilename = newimage->Dir()+"/kernel_from_"+refimage->Name()+".dat";
    //  obj_output<xmlstream> oo(outputfilename);
    double phoratio = 1;
    double chi2 = 0.;
    int nstars = 0;
    int nparams = 0;
    
    
    
    if(makesub) {
      cout << "image_sub : creating sub " << refimage->Name() << " " << newimage->Name() << endl;
      ImageSubtraction sub("sub",refimage,newimage, /*noswap = */  true );
      sub.MakeFits();
      cout << "image_sub : writing kernel in " << outputfilename << " ..." << endl; 
      const KernelFit& kernelFit = sub;
      kernelFit.write(outputfilename);

      phoratio = sub.KernAtCenterSum();
      chi2 = sub.Chi2();
      nstars = sub.NStampsUsed();
      nparams = sub.NParams();
      
    }else{
      cout << "image_sub : psfmatch "  << refimage->Name() << " " << newimage->Name() << endl;
      KernelFitter psfmatch(refimage,newimage,true);
      if( ! psfmatch.DoTheFit()){
	cout << "psfmatch failed" << endl;
	
	// dump something for log
	return EXIT_FAILURE;
      }
      cout << "image_sub : writing kernel in " << outputfilename << " ..." << endl; 
      //    obj_output<xmlstream> oo(outputfilename);
      psfmatch.write(outputfilename);
      //    oo << *(psfmatch.GetKernelFit());
      
      phoratio = psfmatch.PhotomRatio();
      chi2 = psfmatch.Chi2();
      nstars = psfmatch.NStampsUsed();
      nparams = psfmatch.NParams();
    }
    
    FitsHeader head(newimage->FitsName(), RW);
    head.AddOrModKey("PMRATIO", phoratio);
    head.AddOrModKey("PMCHI2" , chi2);
    head.AddOrModKey("PMNSTAR", nstars);
    head.AddOrModKey("PMNPAR" , nparams);
    
    cout << "the end" << endl;
  }catch(PolokaException p)
    {
      p.PrintMessage(cout);
      return EXIT_FAILURE;	      
    } 
  return EXIT_SUCCESS;
}
