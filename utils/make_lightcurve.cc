#include "lightcurveguru.h"

static void usage(const char *pgname)
{
  cerr << pgname << " <filename> " << endl ;
}

int main(int argc, char **argv)
{
  if (argc < 2)  {usage(argv[0]);  exit(1);}

  string file = argv[1];

  LightCurveGuru guru;

  if (!guru.read_init(file))
    {
      cerr << argv[0] << " : Error reading file " << file << endl;
      return EXIT_FAILURE;
    }  

  guru.MonopolizeCPU();

  return EXIT_SUCCESS;
}

