#include "lightcurveguru.h"



/*! Dump help for this program */
void usage(const char* programname) {
  
  cout << " syntax  : " << programname << " <dbimage> " << endl;
  cout << " options : " << endl;
  cout << "            --catalogue (-c) # : path to a catalogue of standard stars (fiducials) " << endl;
  exit(-1);
}


int main(int argc, char **argv)
{
  if (argc < 2) {
    usage(argv[0]);
  }
  bool usecata = false;
  string cata = "";
  string config_file = "";
  
  for (int i=1; i < argc; i++) {
    if (!strcmp(argv[i],"--catalogue") || !strcmp(argv[i],"-c")) {
      usecata = true;
      cata = argv[++i];
      continue;
    }
    if (!strcmp(argv[i],"--help") || !strcmp(argv[i],"-h")) {
      usage(argv[0]);
    }
    config_file = argv[i];
  }
  LightCurveGuru guru;
  if (!guru.read(config_file))
    {
      cerr << " FAILURE reading the lightcurve file " << endl;
      return EXIT_FAILURE;
    }  
  if(usecata)
    guru.UseStdCatalogue(cata);
  
  if (!guru.MonopolizeCPU()) return EXIT_FAILURE;
  cout << " All done " << endl;
  return EXIT_SUCCESS;
}

