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
  
  Kernel galaxy;
  string galname = lc.Ref->name+"/galaxy_sn.fits";
  if(!FileExists(galname)) {
    cout << " cannot open " << galname << endl;
    return -12;
  }
  
  galaxy.readFits(galname);
  int hx = galaxy.HSizeX();
  int hy = galaxy.HSizeY();
  

  // get sn position and max flux
  float x,y,fluxmax;
  {
    string lcname = lc.Ref->name+"/lightcurve_sn.dat";
    if(!FileExists(lcname)) {
      cout << " cannot open " << lcname << endl;
      return -12;
    }
    ifstream lcfile(lcname.c_str());
    string line;
    float flux;
    fluxmax = -12;
    while(getline(lcfile,line)) {
      if(line[0]=='#')
	continue;
      sscanf(line.c_str(),"x %f y %f flux %f",&x,&y,&flux);
      if(flux>fluxmax)
	fluxmax=flux;
    }
    lcfile.close();
  }
  
  cout << "x y fluxmax = " << x << " " << y << " " << fluxmax << endl;
  
  // update PSF
  vignet->Resize(hx,hy);
  vignet->Psf.Tabulate(Point(x,y),*(vignet->psf),(const Window&)*vignet);
  
  
  cout << "hx hy = " << hx << " " << hy << endl;
  
  // write image
  Kernel model = galaxy;
  for(int j=-hy;j<=hy;++j) {
      for(int i=-hx;i<=hx;++i) {
	model(i,j)+=fluxmax*vignet->Psf(i,j);
      }
  }
  string modelname = lc.Ref->name+"/model_sn.fits";
  model.writeFits(modelname);
  

  return EXIT_SUCCESS;
}

