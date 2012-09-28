#include <iostream>
#include "fitsimage.h"
#include "apersestar.h"
#include "reducedimage.h"
#include "allreducedimage.h"
#include "fastfinder.h"
#include "frame.h"
#include "imageutils.h"
#include "vutils.h"
#include "polokaexception.h"

using namespace std;

void usage(const char* progName) {
  cerr << progName << " -c CALIBRATIONFILE DBIMAGE...\n";
  exit(EXIT_FAILURE);
}

static double sqr(const  double& x) { return x*x; }

// compute median and M.A.D. = median(|x - median(x)|)
// robust estimator of standard deviation
static double median_mad(vector<double>& x, double& disp) {
  size_t n = x.size();
  sort(x.begin(), x.end());
  double med = (n & 1) ? x[n/2] : (x[n/2-1] + x[n/2])*0.5;  
  for (vector<double>::iterator it = x.begin(); it != x.end(); ++it) {
    *it = fabs(*it - med);
  }
  sort(x.begin(), x.end());
  double mad = (n & 1) ? x[n/2] : (x[n/2-1] + x[n/2])*0.5;  
  disp = 1.4826 * mad;
  return med;
}

struct ImageZeroPoint {
  
  BaseStarList calibStars;

  void operator () (ReducedImageRef Im) const {

    GtransfoRef raDec2Pix = Im->RaDecToPixels();
    GtransfoRef pix2RaDec = Im->PixelsToRaDec();
    BaseStarList calStars;
    calibStars.ExtractInFrame(calStars, ApplyTransfo(Im->UsablePart(), *pix2RaDec).Rescale(1.1));

    //PSFStarList measStars(Im->Dir()+"/psfstars.list");
    AperSEStarList measStars(Im->AperCatalogName());
    vector<double> zparray;
    FastFinder finder(*(BaseStarList*)&measStars);

    //int nzp = 0;
    for (BaseStarCIterator it = calStars.begin(); it!= calStars.end(); ++it) {
      const BaseStar *calstar = *it;
      const BaseStar *measstar = finder.FindClosest(raDec2Pix->apply(*calstar), 2);
      if (measstar && measstar->flux > 0 && measstar->eflux > 0) {
	zparray.push_back(calstar->flux + 2.5*log10(measstar->flux));
	//++nzp;
      }
    }
    
    double zprms;
    //double zp = clipmean(zparray, nzp, zprms, 3, 5);
    //delete [] zparray;
    double zp = median_mad(zparray, zprms);
    double oldzp = Im->AnyZeroPoint();
    cout << " ImagePSFZeroPoint: " << Im->Name() << " old zp = " << oldzp << endl;
    cout << " ImagePSFZeroPoint: " << Im->Name() << " new zp = " << zp << " +/- " << zprms
	 << " nstars " << zparray.size() << endl;
    FitsHeader head(Im->FitsName(), RW);
    head.AddOrModKey("ZP_PHOT", zp, "Zero Point from measured stars flux matched to catalog");
    head.AddOrModKey("EZP_PHOT", zprms, "R.M.S of ZP_PHOT");
    head.AddOrModKey("NZP_PHOT", int(zparray.size()), "Number of stars used to compute ZP_PHOT");
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
