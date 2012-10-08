#include "fitsimage.h"

static bool fitsbits_process(int bitpix, const char* fitsname) {
  
  FitsImage image(fitsname);
  if (!image.IsValid()) {
    cerr << " invalid FITS file : "  << fitsname << endl;
    return false;
  }
  int oldbitpix = image.KeyVal("BITPIX");
  if (bitpix == oldbitpix) return true;
  // we have to remove and rewrite the file to make sure the changed bitpix 
  // is actually applied. (try fitsbits 16 test.fits on a -32 image 
  // and compare initial and new file size)
  // not sure whether bug is in cfitsio or fitsimage.cc
  remove(fitsname);
  FitsImage newimage(fitsname, image, image);
  newimage.ModKey("BITPIX", bitpix);
  return true;
}

static void usage(const char *progname) {
  cerr << "Usage: " << progname << " BITPIX FITS...\n"
       << "Convert a FITS image with a new BITPIX\n";
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {

  if(argc <=1) usage(argv[0]);
  
  bool ok = true;
  int bitpix = atoi(argv[1]);

  for (int i=2; i< argc; ++i)
    ok &= fitsbits_process(bitpix, argv[i]);
  
  return (ok ? EXIT_SUCCESS: EXIT_FAILURE);
}
