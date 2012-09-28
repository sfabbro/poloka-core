#include <iostream>
#include <fstream>
#include "fitsimage.h"
#include "psfstar.h"
#include "reducedimage.h"
#include "allreducedimage.h"
#include "fastfinder.h"
#include "daophotio.h"
#include "frame.h"
#include "imageutils.h"
#include "vutils.h"
#include "polokaexception.h"

using namespace std;

void usage(const char* progName) {
  cerr << progName << " DBIMAGE...\n";
  exit(EXIT_FAILURE);
}

static double sqr(const  double& x) { return x*x; }

struct ImageZeroPoint {

  void operator () (ReducedImageRef Im) const {
    FitsHeader head(Im->FitsName(),RW);
    double zp =
      double(head.KeyVal("PHOT_C")) + 
      double(head.KeyVal("PHOT_K")) * (double(head.KeyVal("AIRMASS")) - 1) +
      2.5 * log10(double(head.KeyVal("EXPTIME")));
    cout << Im->Name() << ": old zp " << head.KeyVal("ZP_PHOT") << " new zp " << zp << endl;
    head.AddOrModKey("ZP_PHOT", zp, "Zero Point from Elixir PHOT_C see formula");
  }
};

  
int main(int nargs, char **args) {

  if (nargs < 2) usage(args[0]);
  
  try {
  
    ReducedImageList imList;
    ImageZeroPoint imZeroPoint;

    for (int i=1; i<nargs; ++i) {
      char *arg = args[i];
      if (arg[0] != '-') {
	ReducedImageRef im = ReducedImageNew(arg);
	if (!im || !im->IsValid()) { 
	  cerr << args[0] << ": " << arg << " is not a valid dbimage\n";
	  continue;
	}
	imList.push_back(im);
	continue;
      }

      switch (arg[1]) {
      default: usage(args[0]);
      }
    }
    
    if (imList.empty()) throw("missing a valid dbimage");
    
    for_each(imList.begin(), imList.end(), imZeroPoint);

  }  catch (PolokaException(e)) {
    e.PrintMessage(cerr);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
