#include <iostream>
#include <fstream>
#include "fitsimage.h"
#include "reducedimage.h"
#include "allreducedimage.h"
#include "fastfinder.h"
#include "sestar.h"
#include "daophotio.h"
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

struct CompareZeroPoint {
  
  BaseStarList calibStars;

  void operator () (ReducedImageRef Im) const {
    FitsHeader head(Im->FitsName());
    double zp =
      2.5 * log10(double(head.KeyVal("EXPTIME")))
      + double(head.KeyVal("PHOT_C"))
      + double(head.KeyVal("PHOT_K")) * (double(head.KeyVal("AIRMASS")) - 1);

    GtransfoRef raDec2Pix = Im->RaDecToPixels();
    GtransfoRef pix2RaDec = Im->PixelsToRaDec();
    BaseStarList calStars;
    calibStars.ExtractInFrame(calStars, 
			      ApplyTransfo(Im->UsablePart(),
					   *pix2RaDec).Rescale(1.1));

    SEStarList sexStars(Im->CatalogName());    
    FastFinder sfinder(*(BaseStarList*)&sexStars);

    ofstream msex((Im->Dir()+"/zpsex.list").c_str());
    vector<double> dmags;
    for (BaseStarCIterator it = calStars.begin(); it!= calStars.end(); ++it) {
      const BaseStar *cstar = *it;
      const BaseStar *sstar = sfinder.FindClosest(raDec2Pix->apply(*cstar), 0.5);
      if (sstar && sstar->flux > 0 && sstar->eflux > 0) {
	msex << cstar->flux << " "
	     << -2.5*log10(sstar->flux) + zp << endl;
	dmags.push_back(cstar->flux + 2.5*log10(sstar->flux) - zp);
      }
    }
    double zpe;
    cout << Im->Name() << " " 
	 << median_mad(dmags, zpe);
    
    if (!FileExists(Im->Dir()+"/calibrated.als")) { 
      cout << " " << zpe << endl;
      return;
    }

    DaoStarList daoStars;
    ReadDaoList(Im->Dir()+"/calibrated.als", daoStars);
    FastFinder dfinder(*(BaseStarList*)&daoStars);

    ofstream mout((Im->Dir()+"/zpmatch.list").c_str());
    ofstream mlst((Im->Dir()+"/zpdao.lst").c_str());
    write_dao_header<AllstarAls>(mlst, *Im);
    vector<double> dmagd, dmagc;
    for (BaseStarCIterator it = calStars.begin(); it!= calStars.end(); ++it) {
      const BaseStar *cstar = *it;
      const BaseStar *sstar = sfinder.FindClosest(raDec2Pix->apply(*cstar), 0.5);
      const BaseStar *dstar = dfinder.FindClosest(raDec2Pix->apply(*cstar), 0.5);
      if (sstar && sstar->flux > 0 && sstar->eflux > 0) {
	mout << cstar->flux << " "
	     << -2.5*log10(sstar->flux) + zp << " ";
	if (dstar && dstar->flux > 0 && dstar->eflux > 0) {
	  mout << -2.5*log10(dstar->flux) + zp << " ";
	  dmagd.push_back(cstar->flux + 2.5*log10(dstar->flux) - zp);
	  dmagc.push_back(-2.5*log10(dstar->flux/sstar->flux));
	  write_dao_star<AllstarAls>(mlst,
				     *(dynamic_cast<const DaoStar*>(dstar)));
	  mlst << endl;
	} else
	  mout << " 0 0 ";
	mout << endl;
      }
    }
    double mad;
    cout << " " << median_mad(dmagd, mad)
	 << " " << median_mad(dmagc, mad)
	 << " " << zpe
	 << endl;
  }
};

  
int main(int nargs, char **args) {

  if (nargs < 2) usage(args[0]);
  
  try {
  
    ReducedImageList imList;
    CompareZeroPoint cmpZeroPoint;

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
      case 'c': { cmpZeroPoint.calibStars.read(args[++i]); continue; }
      default: usage(args[0]);
      }
    }
    
    if (imList.empty()) 
      throw("missing a valid dbimage");
    if (cmpZeroPoint.calibStars.empty()) 
      throw("calibration catalog is empty");
    
    for_each(imList.begin(), imList.end(), cmpZeroPoint);

  } catch (PolokaException(e)) {
    e.PrintMessage(cerr);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
