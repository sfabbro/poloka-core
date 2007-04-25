#include <iostream>
#include <fstream>
#include <string>
#include <fileutils.h>
#include <fitsimage.h>
#include <gtransfo.h>
#include <wcsutils.h>


using namespace std;

static void usage(const char *pgname)
{
  cerr << pgname << " <vignet_fits> <ref_fits> <x_of_vignet_center> <y_of_vignet_center>" << endl ;
}

int main(int argc, char **argv)
{
  if (argc != 5)  {usage(argv[0]);  exit(1);}
  
  string filename;
  filename = argv[1]; if(! FileExists(filename)) { cerr << "FATAL " << filename << " does not exists" << endl; return EXIT_FAILURE;}
  FitsHeader vignet(filename,RW);
  filename = argv[2]; if(! FileExists(filename)) { cerr << "FATAL " << filename << " does not exists" << endl; return EXIT_FAILURE;}
  FitsHeader ref(filename);
  
  float x  = atof(argv[3]);
  float y  = atof(argv[4]);
  int hsize = (int(vignet.KeyVal("NAXIS1"))-1)/2;
  //int hsize = (vignet.Nx()-1)/2;
  
  TanPix2RaDec TanWcs;
  if(! TanWCSFromHeader(ref,TanWcs)) {
    cerr << "TanWCSFromHeader failed" << endl;
    return EXIT_FAILURE;
  }
  
  GtransfoLinShift shift(x-hsize,y-hsize);
  TanWCS2Header(vignet,TanWcs*shift);
  
  

  return EXIT_SUCCESS;
}

