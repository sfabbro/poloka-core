#include <iostream>
#include "fitsimage.h"
#include "psfstar.h"
#include "reducedimage.h"
#include "allreducedimage.h"
#include "fastfinder.h"
#include "frame.h"
#include "imageutils.h"
#include "vutils.h"
#include "polokaexception.h"

using namespace std;

void usage(const char* progName) {
  cerr << progName << "-c CALIBRATIONFILE DBIMAGE...\n";
  exit(EXIT_FAILURE);
}

static double sqr(const  double& x) { return x*x; }

struct ImageZeroPoint {
  
  BaseStarList calibStars;

  void operator () (ReducedImageRef Im) const {

    GtransfoRef raDec2Pix = Im->RaDecToPixels();
    GtransfoRef pix2RaDec = Im->PixelsToRaDec();
    BaseStarList calStars;
    calibStars.ExtractInFrame(calStars, ApplyTransfo(Im->UsablePart(), *pix2RaDec).Rescale(1.1));

    PSFStarList psfStars(Im->Dir()+"/psfstars.list");
    double *zparray = new double[psfStars.size()];
    double *pzp = &zparray[0];
    FastFinder finder(*(BaseStarList*)&psfStars);

    int nzp = 0;
    for (BaseStarCIterator it = calStars.begin(); it!= calStars.end(); ++it) {
      const BaseStar *calstar = *it;
      const BaseStar *psfstar = finder.FindClosest(raDec2Pix->apply(*calstar), 2);

      if (psfstar && psfstar->flux > 0 && psfstar->eflux > 0) {
	*pzp++ = calstar->flux + 2.5*log10(psfstar->flux);
	++nzp;
      }
    }
    
    double zprms;
    double zp = clipmean(zparray, nzp, zprms, 3, 5);
    delete [] zparray;
    cout << " ImagePSFZeroPoint: " << Im->Name() << " old zp = " << Im->ZP() << endl;
    cout << " ImagePSFZeroPoint: " << Im->Name() << " new zp = " << zp << " +/- " << zprms
	 << " nstars " << nzp << endl;
    FitsHeader head(Im->FitsName(),RW);
    head.AddOrModKey("ZP_PHOT", zp, "Zero Point from PSF stars matched to catalog");
    head.AddOrModKey("EZP_PHOT", zprms, "R.M.S of ZP_PHOT");
    head.AddOrModKey("NZP_PHOT", nzp, "Number of stars used to compute ZP_PHOT");
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
      case 'c': { imZeroPoint.calibStars.read(args[++i]); continue; }
      default: usage(args[0]);
      }
    }
    
    if (imList.empty()) throw("missing a valid dbimage");
    if (imZeroPoint.calibStars.empty()) ("calibration catalog is empty");
    
    for_each(imList.begin(), imList.end(), imZeroPoint);

  }  catch (PolokaException(e)) {
    e.PrintMessage(cerr);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
