#include <string>
#include <iostream>
#include <fstream>

#include "allreducedimage.h"
#include "swarpstack.h"
#include "fitsimage.h"

static void usage(const char* ProgName)
{
  cerr << ProgName << " [OPTION]...FILE\n"
       << "Shift and stack with SWarp\n"
       << "     -o DBIMAGE : output stack DBIMAGE (default: shiftswarp)\n"
       << "     -c FILE    : specify a SWarp configuration file instead of default\n\n";
  exit(EXIT_FAILURE);
}

static void wcs_shift(const string& fitsname, const double& dx, const double& dy) {
  FitsHeader head(fitsname, RW);
  double crpix1 = head.KeyVal("CRPIX1");
  double crpix2 = head.KeyVal("CRPIX2");
  head.AddOrModKey("CRPIX1",crpix1 - dx, "Changed for shift & stack");
  head.AddOrModKey("CRPIX2",crpix2 - dy, "Changed for shift & stack");
}

static void wcs_unshift(const string& fitsname, const double& dx, const double& dy) {
  FitsHeader head(fitsname, RW);
  double crpix1 = head.KeyVal("CRPIX1");
  double crpix2 = head.KeyVal("CRPIX2");
  head.AddOrModKey("CRPIX1",crpix1 + dx,"Original value before shift and stack");
  head.AddOrModKey("CRPIX2",crpix2 + dy,"Original value before shift and stack");
}

int main(int nargs, char **args) {

  if (nargs < 2) usage(args[0]);
  string outName = "shiftswarp";
  string cardsName;
  const char* fileName;

  for (int i=1; i<nargs; ++i) {
    char *arg = args[i];
    if (arg[0] == '-') {
      switch (arg[1]) {
      case 'o': outName = args[++i]; break;
      case 'c': cardsName = args[++i]; break;
      default: usage(args[0]);
      }
      continue;
    }
    fileName = arg;
  }
    
  ReducedImageList ril;
  ifstream ifs(fileName);
  if (ifs.is_open()) {
    char c;
    while (ifs >> c) {
      ifs.unget();
      string name;  double mjd, dx, dy;
      ifs >> name >> mjd >> dx >> dy;
      ReducedImageRef ri = ReducedImageNew(name);
      if (ri && ri->IsValid()) {
	wcs_shift(ri->FitsName(), dx, dy);
	wcs_shift(ri->FitsWeightName(), dx, dy);
	ril.push_back(ri);
      } else 
	cerr << " not a valid dbimage: " << name << endl;
    }
    ifs.close();
  }
  
  if (ril.empty()) {
    cerr << " no input images provided\n";
    return EXIT_FAILURE;
  }

  if (outName.empty()) {
    cerr << " no output name provided\n";
    return EXIT_FAILURE;
  }
  
  if (!cardsName.empty()) SetSwarpCardsName(cardsName);
  
  ReducedImageRef ref;  Frame frame;
  SwarpStack ss(outName, ril, ref, frame);
  ss.MakeFits();
  ss.MakeSatur();
  ss.MakeCatalog();

  ifs.open(fileName);
  if (ifs.is_open()) {
    char c;
    while (ifs >> c) {
      ifs.unget();
      string name;  double mjd, dx, dy;
      ifs >> name >> mjd >> dx >> dy;
      ReducedImage ri(name);
      if (ri.IsValid()) {
	wcs_unshift(ri.FitsName(), dx, dy);
	wcs_unshift(ri.FitsWeightName(), dx, dy);
      } else 
	cerr << " not a valid dbimage: " << name << endl;
    }
    ifs.close();
  }

  return EXIT_SUCCESS;
}
