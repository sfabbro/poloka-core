#include <algorithm>

#include "allreducedimage.h"
#include "imagesubtraction.h"
#include "polokaexception.h"

static void usage(const char *progName) {
  cerr << "Usage: " << progName << " [OPTION]...[DBIMAGE]...\n"
       << "PSF match and subtract DBIMAGE to first DBIMAGE\n\n"
       << "   -d: perform candidate detection on each subtraction\n"
       << "   -f: force first DBIMAGE to be convolved even if worse seeing\n"
       << "   -n: do not subtract, only PSF match\n"
       << "   -o: overwrite\n";
  exit(EXIT_FAILURE);
}

struct ImageSubtract {

  bool overwrite, noswap, dosub, dodetect;
  ReducedImageRef Ref;

  ImageSubtract() : overwrite(false), noswap(false), dosub(true), dodetect(false) {}
  
  void operator () (const ReducedImageRef Im) const {

    // read or redo kernel fit
    KernelFitter kernfit(Ref, Im, noswap);
    if (overwrite || !kernfit.ReadKernel()) {
      if (kernfit.DoTheFit())
	kernfit.WriteKernel(overwrite);
      else {
	PolokaException("PSF match between " + 
			Ref->Name() + " and " + 
			Im->Name() + " failed");
	return;
      }
    }

    if (!dosub) return;
    ImageSubtraction sub(SubtractedName(Ref->Name(), Im->Name()), Ref, Im, kernfit);

    if (overwrite) { 
      if (sub.HasImage()) remove(sub.FitsName().c_str());
      if (sub.HasWeight()) remove(sub.FitsWeightName().c_str());
      if (sub.HasCosmic()) remove(sub.FitsCosmicName().c_str());
      if (FileExists(sub.DetectionsName())) remove(sub.DetectionsName().c_str());
    }

    sub.Execute(DoFits|DoWeight|DoCosmic);
    if (dodetect) sub.MakeCatalog();
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
      im->Execute(DoFits|DoCatalog|DoAperCatalog);
      imList.push_back(im);
      continue;
    }
    switch (arg[1]) {
    case 'd': imSubtract.dodetect = true; break;
    case 'f': imSubtract.noswap = true; break;
    case 'n': imSubtract.dosub = false; break;
    case 'o': imSubtract.overwrite = true; break;
    default : usage(args[0]);
    }
  }

  imSubtract.Ref = imList.front();
  imList.pop_front();
  for_each(imList.begin(), imList.end(), imSubtract);

  return EXIT_SUCCESS;
}
