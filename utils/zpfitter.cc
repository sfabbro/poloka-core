#include <iostream>
#include <fstream>
#include <fileutils.h>
#include <dictfile.h>
#include <basestar.h>
#include <vutils.h>

#include <map>
#include <iomanip>

using namespace std;

static void usage(const char *pgname)
{
  cerr << pgname << " -c <dictfile> -m <keyofreferencemagnitude> -f <keyofflux>" << endl ;
  exit(1);
}

int main(int argc, char **argv)
{
  
  string catalogname = "";
  string magkey = "";
  string fluxkey = "";
  
  if (argc != 7)  {usage(argv[0]);}
  for (int i=1; i<argc; ++i)
    {
      char *arg = argv[i];
      if (arg[0] != '-')
	{
	  cerr << "unexpected parameter " << arg << endl;
	  usage(argv[0]);
	}
      switch (arg[1])
	{
	case 'm' : magkey = argv[++i]; break;
	case 'c' : catalogname = argv[++i]; break;
	case 'f' : fluxkey = argv[++i]; break;
	default : 
	  cerr << "unknown option " << arg << endl;
	  usage(argv[0]);
	}
    }
  
  
  if(!FileExists(catalogname)) {
    cerr << "cant find catalog " << catalogname << endl;
    usage(argv[0]);
  }
  
  DictFile catalog(catalogname);
  if ( ! catalog.HasKey(magkey) ) {
    cerr << "catalog " << catalogname << " does have any key '" << magkey << "'" << endl;
    usage(argv[0]);
  }
  if ( ! catalog.HasKey(fluxkey) ) {
    cerr << "catalog " << catalogname << " does have any key '" << fluxkey << "'" << endl;
    usage(argv[0]);
  }
  
  // fill a vector of ZPs
  unsigned int nentries = catalog.size();
  double* values = new double[nentries];
  double flux,mag;
  int count = 0;
  for(DictFileCIterator entry=catalog.begin();entry!=catalog.end();++entry) {
    
    flux = entry->Value(fluxkey);
    mag = entry->Value(magkey);
    
    if(flux<=0) continue;
    values[count++]=2.5*log10(flux)+mag;
  }
  if (count==0) {
    cerr << "no valid entries in catalog!" << endl;
    delete [] values;
    return EXIT_FAILURE;
  }
  double zp,rms;
  FILE *file = fopen("zpfitter.dat","w");
  
  // gaussian
  zp = gaussianfit(values,count,zp,rms,3.,true);
  zp = gaussianfit(values,count,zp,rms,2.,false);
  zp = gaussianfit(values,count,zp,rms,1.5,false);
  zp = gaussianfit(values,count,zp,rms,1.,false);
  fprintf(file,"gaussian_zp_rms_count= %6.6f %6.6f %d\n",zp,rms,count);
  printf("gaussian_zp_rms_count= %6.6f %6.6f %d\n",zp,rms,count);
  fclose(file);
  delete [] values;
  return EXIT_SUCCESS;
}


