#include <iostream>
#include <fstream>
#include <string>
#include <lightcurve.h>
#include <simfitphot.h>
#include <dimage.h>

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
  LightCurve &lc = *(fids.begin());
  int radius = 20;
  
  PhotStar* star = lc.Ref;
  SimFitRefVignet *vignet = new SimFitRefVignet(star,lc.Ref->Image(),radius);
 
  int hx,hy;
  {
    // get size of reference data vignet 
    Kernel refdata;
    string name = vignet->Image()->Name();
    string refdataname = "T_"+name+name.substr(0,9)+"_sn_data.fits";
    if(!FileExists(refdataname)) {
      cout << " cannot open " << refdataname << endl;
      return -12;
    }
    refdata.readFits(refdataname);
    hx = refdata.HSizeX();
    hy = refdata.HSizeY();
    cout << "size used for vignet hx hy = " << hx << " " << hy << endl;
  }
  

  // get sn position
  float x,y;
  {
    string lcname = "lightcurve_sn.dat";
    if(!FileExists(lcname)) {
      cout << " cannot open " << lcname << endl;
      return EXIT_FAILURE;
    }
    ifstream lcfile(lcname.c_str());
    string line;
    float flux;
    while(getline(lcfile,line)) {
      if(line[0]=='#')
	continue;
      if(sscanf(line.c_str(),"%f %f %f",&x,&y,&flux)!=3)
	if(sscanf(line.c_str(),"x %f y %f flux %f",&x,&y,&flux)!=3) {
	  cerr << " problem when reading " << lcname << endl;
	  cerr << " must quit " << endl;
	  return EXIT_FAILURE;
	}
    }
    lcfile.close();
  }

  // get max flux
  float fluxmax=0;
  if(true){
    string lcname = "lc2fit_per_night.dat";
    if(!FileExists(lcname)) {
      cout << " cannot open " << lcname << endl;
      return EXIT_FAILURE;
    }
    ifstream lcfile(lcname.c_str());
    string line;
    fluxmax = -12;
    float date,flux,fluxerror,zp;
    while(getline(lcfile,line)) {
      if(line[0]=='#')
	continue;
      if(line[0]=='@')
	continue;
      if(sscanf(line.c_str(),"%f %f %f %f",&date,&flux,&fluxerror,&zp)==4) {
	if(flux>fluxmax)
	  fluxmax=flux;
      }
      
    }
    lcfile.close();
  }
  
  cout << "x y fluxmax = " << x << " " << y << " " << fluxmax << endl;
  
  // read galaxy
  {
    Kernel galaxy;
    string galname = "galaxy_sn.fits";
    if(!FileExists(galname)) {
      cout << " cannot open " << galname << endl;
      return -12;
    }
    galaxy.readFits(galname);
    Kernel smallgal(hx,hy);
    for(int j=-hy;j<=hy;++j) {
      for(int i=-hx;i<=hx;++i) {
	smallgal(i,j)=galaxy(i,j);
      }
    }
    smallgal.writeFits("model_gal.fits");
  }
  
  // update PSF
  vignet->Resize(hx,hy);
  vignet->Psf.Tabulate(Point(x,y),*(vignet->psf),(const Window&)*vignet);
  for(int j=-hy;j<=hy;++j) {
      for(int i=-hx;i<=hx;++i) {
	vignet->Psf(i,j)*=fluxmax;
      }
    }
  vignet->Psf.writeFits("model_psf.fits");
  
  

  return EXIT_SUCCESS;
}

