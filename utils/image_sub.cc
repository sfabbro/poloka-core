#include <iostream>
#include <vector>

#include "imagesubtraction.h"
#include "fitsimage.h"

// for io
#include "kernelfit_dict.h"
#include "objio.h"
#include "typemgr.h"


static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " ref new  " << endl;
  cerr << "  match PSF of <DbImage(s)> with a reference image \n" 
       << "   Note:both images have to be on the same pixel grid ! \n";

  exit(-1);
}

int main(int nargs, char **args)
{
  // if nothing is given
  if (nargs < 3){usage(args[0]);}
  
  ReducedImage refimage(args[1]);
  ReducedImage newimage(args[2]);
  string subname = "newsub3";
  ImageSubtraction sub(subname,refimage,newimage);
  string outputfilename = "kernel_from_"+sub.Best()->Name()+"_to_"+sub.Worst()->Name()+".xml";
  bool kernelissaved = FileExists(outputfilename);
  KernelFit *kernel=0;
  if(kernelissaved) {
    cout << "image_sub : reading kernel from " << outputfilename << " ..." << endl;
    kernel = new KernelFit();
    obj_input<xmlstream> oi(outputfilename);
    oi >> *kernel;
    oi.close();
    sub.SetKernelFit(kernel);
    cout << "-----------------------------------" << endl;
    cout << *(sub.GetKernelFit()) << endl;
    cout << "-----------------------------------" << endl;
    //return -1;
  }
  sub.Execute(DoFits| DoWeight | DoCatalog);
  if(!kernelissaved) {
    cout << "image_sub : writing kernel in " << outputfilename << " ..." << endl;
    obj_output<xmlstream> oo(outputfilename);
    oo << *(sub.GetKernelFit());
    oo.close();
    cout << "-----------------------------------" << endl;
    cout << *(sub.GetKernelFit()) << endl;
    cout << "-----------------------------------" << endl;
  }
  cout << "the end" << endl;

  return 0;
}
