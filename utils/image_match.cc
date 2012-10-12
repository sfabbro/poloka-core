#include <string>

#include "reducedutils.h"
#include "allreducedimage.h"
#include "subimage.h"
#include "transformedimage.h"
#include "polokaexception.h"

static void usage(const char *progName) {
  cerr << "Usage: " << progName << " [OPTION]...DBIMAGE...\n"
       << "Match DBIMAGE to a geometric reference image (first DBIMAGE)\n\n" 
       << "   -i : integer shifting (no interpolation)\n"
       << "   -n : no resampling, only match catalogues\n"
       << "   -t x y: translation parameters\n"
       << "   -u : union of all frames instead of intersection\n"
       << "   -w : match using WCS information only\n\n";
  exit(EXIT_FAILURE);
}

struct ImageMatcher {

  ImageMatcher() : doResample(true), doIntShift(false), wcsOnly(false) {}

  bool doResample, doIntShift, wcsOnly;
  ReducedImageRef Ref;
  GtransfoRef RefToIm, ImToRef;

  void operator () (const ReducedImageRef Im) {
    try {
      if (Ref && *Im == *Ref) {
	cout << " " << Im->Name() << " is same as reference, skipping\n";
	return;
      }
      ReducedImageRef ref = Ref;
      // ugly hack: ref bigger than 3X image, extract subimage
      if (Ref && Ref->XSize()*Ref->YSize() > 3*Im->XSize()*Im->YSize()) {
	cout << " " << Ref->Name() << " is a big image, will extract a sub image using rough match\n";
	ref = new SubImage("Sub_" + Ref->Name() + "_" + Im->Name(), Ref->Name(), Im->Name(), 50);
	ref->Execute(DoFits|DoWeight|DoCatalog);
      }

      if (!ImToRef && !RefToIm && wcsOnly) {
	ImToRef = FindTransfoFromWCS(*Im, *Ref);
	RefToIm = FindTransfoFromWCS(*Ref, *Im);
      }

      if (doResample && ref)
	ImageResample(*Im, *ref, RefToIm, ImToRef);
      else if (doIntShift && ref)
	ImageIntegerShift(*Im, *ref, RefToIm);
      else if (doResample)
	ImageResample(*Im, RefToIm, ImToRef);
      else {
	if (wcsOnly)
	  cout << *ImToRef;
	else
	  cout << *FindTransfo(*Im, *ref);
      }

    } catch (PolokaException p) {
      p.PrintMessage(cerr);
    }
  }
};


int main(int nargs, char **args) {

  if (nargs < 2) usage(args[0]); 

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
    case 'h': usage(args[0]);
    case 'i': imMatcher.doResample = false; imMatcher.doIntShift = true; break;
    case 'n': imMatcher.doResample = false; break;
    case 't': {
      double dx = atof(args[++i]);
      double dy = atof(args[++i]);
      imMatcher.ImToRef = new GtransfoLinShift( dx,  dy);
      imMatcher.RefToIm = new GtransfoLinShift(-dx, -dy);
    }
      break;
    case 'u': doUnion = true; break;
    case 'w': imMatcher.wcsOnly = true; break;
    default : usage(args[0]);
    }
  }

  if (imList.empty()) {
    cerr << args[0] << " error: need at least 1 image";
    return EXIT_FAILURE;
  }

  if (imList.size() == 1 && !imMatcher.ImToRef) {
    cerr << args[0] << " error: need at least a reference and an image";
    return EXIT_FAILURE;
  }

  if (doUnion) {
    string unionRefName = "U_" + imList.front()->Name();
    ReducedImageRef unionRef = ReducedImageNew(unionRefName);
    if (unionRef && unionRef->Execute(ToTransform(*imList.front()))) {
      cout << " Union frame " << unionRefName << " already produced\n";
      imMatcher.Ref = unionRef;
    } else {
      cout << " Creating a union frame reference " << unionRefName << endl;
      MakeUnionRef(imList, *imList.front(), unionRefName, imMatcher.wcsOnly);
      imMatcher.Ref = ReducedImageNew(unionRefName);
    }
    imList.pop_front();
  } else if (!imMatcher.ImToRef && !imMatcher.RefToIm) {
    imMatcher.Ref = imList.front();
    imList.pop_front();
  }

  for_each(imList.begin(), imList.end(), imMatcher);

  return EXIT_SUCCESS;
}
