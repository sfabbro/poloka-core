#include <iostream>
#include <fstream>
#include <fileutils.h>
#include <dicstar.h>
#include <basestar.h>
#include <vutils.h>

#include <map>
#include <iomanip>

using namespace std;

static void usage(const char *pgname)
{
  cerr << pgname << " -c <dicstarlist>" << endl;
  cerr << "options:  -m <keyofreferencemagnitude> (default is \"mag\")"<< endl;
  cerr << "          -f <keyofflux>               (default is \"flux\")" << endl ;
  exit(1);
}

int main(int argc, char **argv)
{
  
  string catalogname = "";
  string magkey = "mag";
  string fluxkey = "flux";
  
  if (argc < 3)  {usage(argv[0]);}
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
  
  DicStarList catalog(catalogname);
  if ( catalog.empty()) {
    cerr << "catalog is empty" << endl;
    usage(argv[0]);
  }
  int nstars=-1;
  if(catalog.GlobVal().HasKey("NSTARS")) nstars=int(catalog.GlobVal().getDoubleValue("NSTARS"));
  int nimages=-1;
  if(catalog.GlobVal().HasKey("NIMAGES")) nimages=int(catalog.GlobVal().getDoubleValue("NIMAGES"));
  
  // fill a vector of ZPs
  unsigned int nentries = catalog.size();
  double* values = new double[nentries];
  double flux,mag;
  int count = 0;
  for(DicStarCIterator entry=catalog.begin();entry!=catalog.end();++entry) {
    
    if(fluxkey=="flux")
      flux = (*entry)->flux;
    else
      flux = (*entry)->getval(fluxkey);
    mag = (*entry)->getval(magkey);
    
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
  fprintf(file,"@PSFZP %6.6f\n",zp);
  fprintf(file,"@PSFZPERROR %6.6f\n",rms);
  fprintf(file,"@NMEASUREMENTS %d\n",count);
  
  printf("@PSFZP %6.6f %6.6f %d\n",zp,rms,count);
  if(nstars>0) {
    fprintf(file,"@NSTARS %d\n",nstars);
  }
  if(nimages>0) {
    fprintf(file,"@NIMAGES %d\n",nimages);
  }
  fclose(file);
  delete [] values;
  return EXIT_SUCCESS;
}


