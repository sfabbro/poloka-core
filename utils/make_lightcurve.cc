#include <iostream>
#include <fstream>
#include <lightcurve.h>
#include <simfitphot.h>

using namespace std;

static void usage(const char *pgname)
{
  cerr << pgname << " <filename> " << endl ;
}

int main(int argc, char **argv)
{
  if (argc < 2)  {usage(argv[0]);  exit(1);}

  ifstream lightfile(argv[1]);
  if (!lightfile) return EXIT_FAILURE;

  LightCurveList fids(lightfile);
  SimFitPhot doFit(fids);
  for_each(fids.begin(), fids.end(), doFit);

  return EXIT_SUCCESS;
}

