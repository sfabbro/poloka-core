#include <iostream>
#include <cmath>

#include <poloka/fitsimage.h>
#include <poloka/sestar.h>
#include <poloka/reducedimage.h>
#include <poloka/fastfinder.h>
#include <poloka/frame.h>
#include <poloka/imageutils.h>
#include <poloka/vutils.h>
#include <poloka/polokaexception.h>


static void usage(const char* progname) {
  cerr << "Usage: " << progName << " -c CALIBRATIONFILE DBIMAGE...\n"
       << "Compute a median zero point for an image, matching to a ra dec mag catalog\n"
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
    calibStars.ExtractInFrame(calStars,
			      ApplyTransfo(Im->UsablePart(),
					   *pix2RaDec).Rescale(1.1));
    
    SEStarList measStars(Im->CatalogName());
    vector<double> zparray;
    FastFinder finder(*(BaseStarList*)&measStars);

    for (BaseStarCIterator it = calStars.begin(); it!= calStars.end(); ++it) {
      const BaseStar *calstar = *it;
      const BaseStar *measstar = finder.FindClosest(raDec2Pix->apply(*calstar), 2);
      if (measstar && measstar->flux > 0 && measstar->eflux > 0) {
	zparray.push_back(calstar->flux + 2.5*log10(measstar->flux));
      }
    }
    
    double zprms;
    double zp = median_mad(zparray, zprms);
    double oldzp = Im->AnyZeroPoint();
    cout << " ImageZeroPoint: " << Im->Name()
	 << " old zp = " << oldzp << endl;
    cout << " ImageZeroPoint: " << Im->Name()
	 << " new zp = " << zp << " +/- " << zprms
	 << " nstars " << zparray.size() << endl;
    FitsHeader head(Im->FitsName(), RW);
    head.AddOrModKey("ZP_PHOT", zp,
		     "Zero Point from measured stars flux matched to catalog");
    head.AddOrModKey("EZP_PHOT", zprms, "R.M.S of ZP_PHOT");
    head.AddOrModKey("NZP_PHOT", int(zparray.size()),
		     "Number of stars used to compute ZP_PHOT");
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
	ReducedImage* im = new ReducedImage(arg);
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
    if (imZeroPoint.calibStars.empty()) throw("calibration catalog is empty");
    
    for_each(imList.begin(), imList.end(), imZeroPoint);

  } catch (PolokaException(e)) {
    e.PrintMessage(cerr);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
