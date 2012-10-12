#include <iostream>

#include "wcsutils.h"
#include "reducedimage.h"
#include "fitsimage.h"

static void usage(const char* ProgName)
{
  cerr << ProgName << " DBIMAGE dx dy\n"
       << "Shift wcs of images\n";
  exit(EXIT_FAILURE);
}

static void wcs_shift(const string& fitsname, const double& dx, const double& dy) {
  FitsHeader head(fitsname, RW);
  GtransfoLinShift tr(dx,dy);
  GtransfoRef wcs = WCSFromHeader(head);
  GtransfoRef wcstr = GtransfoCompose(wcs, &tr);
  TanPix2RaDec *wcstan = dynamic_cast<TanPix2RaDec*>((Gtransfo*)wcstr);
  if (wcstan)
    TanWCS2Header(head, *wcstan);
}

int main(int nargs, char **args) {

  if (nargs < 4) usage(args[0]);

  ReducedImage ri(args[1]);
  double dx = atof(args[2]);
  double dy = atof(args[3]);

  if (ri.IsValid()) {
    wcs_shift(ri.FitsName(), dx, dy);
    wcs_shift(ri.FitsWeightName(), dx, dy);
  } else 
    cerr << " not a valid dbimage: " << args[1] << endl;

  return EXIT_SUCCESS;
}
