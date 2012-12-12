#include <string>
#include <iostream>
#include <fstream>

#include "allreducedimage.h"
#include "swarpstack.h"
#include "fitsimage.h"
#include "astroutils.h"

static void usage(const char* ProgName)
{
  cerr << ProgName << " [OPTION]...FILE\n"
       << "Shift and stack with SWarp\n"
       << "     -o DBIMAGE : output swarped DBIMAGE (default: shiftswarp)\n"
       << "     -c FILE    : specify a SWarp configuration file instead of default\n\n";
  exit(EXIT_FAILURE);
}

static void wcs_shift(const string& fitsname, const double& dcrval1, const double& dcrval2) {
  FitsHeader head(fitsname, RW);
  double crval1 = head.KeyVal("CRVAL1");
  double crval2 = head.KeyVal("CRVAL2");
  head.AddOrModKey("OCRVAL1", crval1, "Original value of CRVAL1");
  head.AddOrModKey("OCRVAL2", crval2, "Original value of CRVAL2");
  head.AddOrModKey("CRVAL1", crval1 + dcrval1, "Changed for shift & stack");
  head.AddOrModKey("CRVAL2", crval2 + dcrval2, "Changed for shift & stack");
}

static void wcs_unshift(const string& fitsname) {
  FitsHeader head(fitsname, RW);
  double crval1 = head.KeyVal("OCRVAL1");
  double crval2 = head.KeyVal("OCRVAL2");
  head.AddOrModKey("CRVAL1", crval1, "Original value");
  head.AddOrModKey("CRVAL2", crval2, "Original value");
  head.RmKey("OCRVAL1");
  head.RmKey("OCRVAL2");
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
      string name;  double dcrval1, dcrval2;
      ifs >> name >> dcrval1 >> dcrval2;
      ReducedImageRef ri = ReducedImageNew(name);
      if (ri && ri->IsValid()) {
	wcs_shift(ri->FitsName(), dcrval1, dcrval2);
	if (ri->HasWeight()) wcs_shift(ri->FitsWeightName(), dcrval1, dcrval2);
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


  for (ReducedImageIterator it=ril.begin(); it!=ril.end(); ++it) {
    ReducedImageRef ri = *it;
    if (ri->IsValid()) {
      wcs_unshift(ri->FitsName());
	if (ri->HasWeight()) wcs_unshift(ri->FitsWeightName());
    } else 
      cerr << " not a valid dbimage: " << ri->Name() << endl;
  }

  return EXIT_SUCCESS;
}
