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
  
  // read galaxy
  Kernel galaxy;
  string galname = "galaxy_sn.fits";
  if(!FileExists(galname)) {
    cout << " cannot open " << galname << endl;
    return -12;
  }
  galaxy.readFits(galname);

  // get size of reference data vignet 
  Kernel refdata;
  string name = vignet->Image()->Name();
  string refdataname = "T_"+name+name.substr(0,9)+"_sn_data.fits";
  if(!FileExists(refdataname)) {
    cout << " cannot open " << refdataname << endl;
    return -12;
  }
  refdata.readFits(refdataname);
  
  // load weights in case of vignet on border of frame
  Kernel weight;
  string weightdataname = "T_"+name+name.substr(0,9)+"_sn_weight.fits";
  if(!FileExists(weightdataname)) {
    cout << " cannot open " << weightdataname << endl;
    return -12;
  }
  weight.readFits(weightdataname);

  
  int hx = refdata.HSizeX();
  int hy = refdata.HSizeY();
  cout << "hx hy = " << hx << " " << hy << endl;


  // get sn position and max flux
  float x,y,fluxmax;
  {
    string lcname = "lightcurve_sn.dat";
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
  
  // create an image with and without sn
  DImage image((2*hx+1)*2+1,(2*hy+1));
  
  int dy;
  float vmax=0;
  float v;
  for(int j=-hy;j<=hy;++j) {
    dy = j+hy;
    for(int i=-hx;i<=hx;++i) {
      if(true) {
	if(weight(i,j)>1.e-10) {
	  image(i+hx,dy)=galaxy(i,j);
	  v = galaxy(i,j) +fluxmax*vignet->Psf(i,j);
	  image(i+hx+(2*hx+1)+1,dy)=v;
	  if(v>vmax)
	    vmax=v;
	}
      }
    }
  }
  for(int j=-hy;j<=hy;++j) {
    image(2*hx+1,hy+j)=vmax;
  }
  
  string modelname = "model_sn.fits";
  image.writeFits(modelname);
  

  return EXIT_SUCCESS;
}

