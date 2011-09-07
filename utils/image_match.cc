#include <string>

#include "reducedutils.h"
#include "allreducedimage.h"
#include "subimage.h"
#include "transformedimage.h"
#include "polokaexception.h"

static void usage(const char *progName) {
  cerr << "Usage: " << progName << " [OPTION]...[DBIMAGE]...\n"
       << "Match DBIMAGE to a geometric reference image (first DBIMAGE)\n\n" 
       << "   -n : no resampling, only match catalogues\n"
       << "   -i : integer shifting (no interpolation)\n"
       << "   -t x y: translation parameters\n"
       << "   -u : union of all frames instead of intersection\n";
  exit(EXIT_FAILURE);
}

struct ImageMatcher {

  ImageMatcher() : doResample(true), doIntShift(false) {}

  bool doResample, doIntShift;
  ReducedImageRef Ref;
  GtransfoRef RefToIm, ImToRef;

  void operator () (const ReducedImageRef Im) const {
    try {
      if (*Im == *Ref) {
	cout << " " << Im->Name() << " is same as reference, skipping\n";
	return;
      }
      ReducedImageRef ref = Ref;
      // ugly hack: ref bigger than 3X image, extract subimage
      if (Ref->XSize()*Ref->YSize() > 3*Im->XSize()*Im->YSize()) {
	cout << " " << Ref->Name() << " is a big image, will extract a sub image using rough match\n";
	ref = new SubImage("Sub_" + Ref->Name() + "_" + Im->Name(), Ref->Name(), Im->Name(), 50);
	ref->Execute(DoFits|DoCatalog);
      }

      if (doResample)
	ImageResample(*Im, *ref, RefToIm, ImToRef);
      else if (doIntShift)
	ImageIntegerShift(*Im, *ref, RefToIm);
      else
	cout << *FindTransfo(*Im, *ref);

    } catch (PolokaException p) {
      p.PrintMessage(cerr);
    }
  }
};


int main(int nargs, char **args) {

  if (nargs < 2) usage(args[0]); 

  if (nargs == 2) {
    cerr << args[0] << " error: need at least 2 images";
    return EXIT_FAILURE;    
  }

  bool doUnion = false;
  ReducedImageList imList;
  ImageMatcher imMatcher;
  
  for (int i=1; i<nargs; ++i) {
    char *arg = args[i];
    if (arg[0] != '-') {
      ReducedImageRef im = ReducedImageNew(arg);
      if (!im || !im->IsValid()) { 
	cerr << " not a valid dbimage: " << arg << endl;
	continue;
      }
      imList.push_back(im);
      continue;
    }
    switch (arg[1]) {
    case 'n': imMatcher.doResample = false; break;
    case 'i': imMatcher.doResample = false; imMatcher.doIntShift = true; break;
    case 't': { 
      double dx = atof(args[++i]);
      double dy = atof(args[++i]);
      imMatcher.ImToRef = new GtransfoLinShift( dx,  dy);
      imMatcher.RefToIm = new GtransfoLinShift(-dx, -dy);
    }
      break;
    case 'u': doUnion = true; break;
    default : usage(args[0]);
    }
  }

  if (doUnion) {
    string unionRefName = "U_" + imList.front()->Name();
    ReducedImageRef unionRef = ReducedImageNew(unionRefName);
    if (unionRef && unionRef->Execute(ToTransform(*imList.front()))) {
      cout << " Union frame " << unionRefName << " already produced\n";
      imMatcher.Ref = unionRef;
    } else {
      cout << " Creating a union frame reference " << unionRefName << endl;
      MakeUnionRef(imList, *imList.front(), unionRefName);
      imMatcher.Ref = ReducedImageNew(unionRefName);
    }
  } else
    imMatcher.Ref = imList.front();

  imList.pop_front();

  for_each(imList.begin(), imList.end(), imMatcher);

  return EXIT_SUCCESS;
}
