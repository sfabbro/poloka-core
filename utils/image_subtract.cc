#include <algorithm>

#include "allreducedimage.h"
#include "imagesubtraction.h"
#include "polokaexception.h"

static void usage(const char *progName) {
  cerr << "Usage: " << progName << " [OPTION]...DBIMAGE...\n"
       << "PSF match and subtract DBIMAGE to first DBIMAGE\n\n"
       << "   -d: perform candidate detection on each subtraction\n"
       << "   -o: output subtraction dbimage name (default: DBIMAGE-FIRST)\n"
       << "   -r: first DBIMAGE will always be convolved even if worse seeing\n"
       << "   -n: do not subtract, only perform and save PSF match\n"
       << "   -f: overwrite\n";
  exit(EXIT_FAILURE);
}

struct ImageSubtract {

  bool overWrite, noSwap, doSub, doDetect;
  ReducedImageRef Ref;
  string subName;

  ImageSubtract() : overWrite(false), noSwap(false), doSub(true), doDetect(false) {}
  
  void operator () (const ReducedImageRef Im) const {

    // read or redo kernel fit
    try {
      KernelFitter kernfit(Ref, Im, noSwap);
      if (overWrite || !kernfit.ReadKernel())
	if (kernfit.DoTheFit())
	  kernfit.WriteKernel(overWrite);
	else
	  cerr << " match between "
	       << Ref->Name() << " and "
	       << Im->Name() << " failed\n";
	  
      if (!doSub) return;
      string subname = subName;
      if (subName.empty()) subname = SubtractedName(Ref->Name(), Im->Name());
      ImageSubtraction sub(subname, Ref, Im, kernfit);
      
      if (overWrite) { 
	if (sub.HasImage()) remove(sub.FitsName().c_str());
	if (sub.HasWeight()) remove(sub.FitsWeightName().c_str());
	if (sub.HasCosmic()) remove(sub.FitsCosmicName().c_str());
	if (FileExists(sub.DetectionsName())) remove(sub.DetectionsName().c_str());
      }
      
      sub.Execute(DoFits|DoWeight|DoCosmic);
      if (doDetect) sub.MakeCatalog();
    } catch (PolokaException p) {
      p.PrintMessage(cerr);
    }
  }
};

int main(int nargs, char **args) {
  if (nargs < 2) usage(args[0]);
  if (nargs == 2) { 
    cerr << args[0] << ": needs an image to subtract from\n";
    return EXIT_FAILURE;
  }

  ReducedImageList imList;
  ImageSubtract imSubtract;
  
  for (int i=1; i<nargs; ++i) {
    char *arg = args[i];
    if (arg[0] != '-') {
      ReducedImageRef im = ReducedImageNew(arg);
      if (!im || !im->IsValid()) { 
	cerr << arg << ": not a valid dbimage\n";
	continue;
      }
      im->Execute(DoFits|DoCatalog);
      imList.push_back(im);
      continue;
    }
    switch (arg[1]) {
    case 'd': imSubtract.doDetect = true; break;
    case 'r': imSubtract.noSwap = true; break;
    case 'n': imSubtract.doSub = false; break;
    case 'f': imSubtract.overWrite = true; break;
    case 'o': imSubtract.subName = args[++i]; break;
    default : usage(args[0]);
    }
  }

  imSubtract.Ref = imList.front();
  imList.pop_front();
  for_each(imList.begin(), imList.end(), imSubtract);

  return EXIT_SUCCESS;
}
