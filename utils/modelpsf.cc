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
      sscanf(line.c_str(),"%f %f %f",&x,&y,&flux);
      if(flux>fluxmax)
	fluxmax=flux;
    }
    lcfile.close();
  }
  
  cout << "x y fluxmax = " << x << " " << y << " " << fluxmax << endl;
  
  // update PSF
  vignet->Resize(hx,hy);
  vignet->Psf.Tabulate(Point(x,y),*(vignet->psf),(const Window&)*vignet);
  vignet->Psf.writeFits("model_psf.fits");
  
  return EXIT_SUCCESS;
}

