#include "vignetphot.h"

static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " [OPTIONS] <DbImage(s)>\n"
       << "  Weighted photometry on detected stars \n"
       << "   OPTIONS are \n"
       << "     -rad VALUE : radius value in FWHM (default is 1 FWHM) \n";

  exit(-1);
}

static void read_phot(PhotStarList &Stars, const string& SexCatName)
{
  SEStarList selist(SexCatName);
  for (SEStarCIterator it=selist.begin(); it != selist.end(); ++it) Stars.push_back(new PhotStar(**it));
}

class processPhot {

  double nfwhm;
  VignetPhot process;

public:

  processPhot(const double& Nfwhm) : nfwhm(Nfwhm) {}

  void operator() (const string& Name)
  {
    ReducedImage rim(Name);
    double fwhm = rim.Seeing()*2.3548;
    GaussPsf psf(rim);
    PhotStarList stars;
    read_phot(stars, rim.CatalogName());
    process.SetRadius(nfwhm*fwhm);
    process.SetSkyRadius(fwhm*2., fwhm*6.);
    process.psf = &psf;
    for_each(stars.begin(), stars.end(), process);
    stars.write("phot.list");
  }

};


int main(int nargs, char **args)
{
  // if nothing is given
  if (nargs < 2) { usage(args[0]); }

  // default stuff
  double rad = 1.;

  StringList imList;

  // loop over arguments
  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      // images
      if (arg[0] != '-') 
	{
	  imList.push_back(arg);
	  continue;
	}

      // options
      sscanf(arg, "%s", ++arg);
      if (strcmp(arg,"rad")==0)   { ++i; rad = atof(args[i]); continue; }

      // unrecognized option
      usage(args[0]);
    }

  processPhot doPhot(rad);

  for_each(imList.begin(), imList.end(), doPhot);

  return EXIT_SUCCESS;
}
