#include <vector>
#include <string>

#include "fileutils.h"
#include "fitsimage.h"
#include "fitstoad.h"
#include "reducedimage.h"
#include "imageutils.h"
#include "wcsutils.h"
#include "astroutils.h"

static bool newhorizon_preprocess(const ReducedImage& Rim) {

  FitsImage *elixir = new FitsImage(Rim.ElixirName());
  if (!elixir->IsValid())
    {
      cerr << " magellan_to_poloka : invalid file : "  << Rim.Name() << endl;
      return false;
    }
  
  //FitsHeader *head = new FitsHeader(Rim.Dir()+"astrogwyn.head");
  //FitsImage image(Rim.FitsName(), *head, *elixir);
  FitsImage image(Rim.FitsName(), *elixir, *elixir);
  delete elixir;
  //delete head;
  return true;
  // trim
  Frame Illu = TotalIlluRegion(image);
  image.Trim(Illu);


  // compute saturation
  Pixel *pim = image.begin();
  Pixel satur = 65535;
  satur = ComputeSaturation(image); 
  image.AddOrModKey("SATURLEV", satur, " Current saturation level");

  // compute sky
  float sky, sig;
  image.SkyLevel(&sky, &sig);
  double rdnoise = image.KeyVal("TOADRDON");
  double gain = image.KeyVal("TOADGAIN");
  double sigth = sqrt(sky*gain + rdnoise*rdnoise) / gain;
  image.AddOrModKey("SKYLEV", sky," Original sky level");
  image.AddOrModKey("SKYSIGEX", sig," Original sky sigma");
  image.AddOrModKey("SKYSIGTH", sigth," app. Theoretical sigma (sqrt(sky*gain+rdnoise^2)/gain");

  // fix RA & DEC
  double ra,dec;
  if (RaDecFromWCS(image, ra, dec)) {
    image.AddOrModKey("RA",  RaDegToString(ra), "Right Ascension from original WCS header");
    image.AddOrModKey("DEC", DecDegToString(dec), "Declination from original WCS header");
  }

  // create saturation image
  FitsHeader &header = image;
  FitsImage sat(Rim.FitsSaturName(), header);
  Pixel *psat = sat.begin();
  for (psat = sat.begin(), pim = image.begin(); pim < image.end() ; ++pim, ++psat) {
    if (*pim >= satur) {
      *psat = 1;
    }
  }
  sat.AddOrModKey("BITPIX",8);

  return true;
}

int main(int argc,char **argv) {
  bool ok = true;
  for (int i=1; i<argc && ok; ++i) {
    ok = newhorizon_preprocess(ReducedImage(argv[i]));
  }
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
