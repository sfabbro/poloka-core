#include <vector>
#include <string>

#include "fileutils.h"
#include "fitsimage.h"
#include "frame.h"
#include "reducedimage.h"
#include "histo1d.h"
#include "superflat.h"
#include "fitstoad.h" 
#include "fringeutils.h"



static void elixir_process(ReducedImage &Rim, bool old_fringe_method, bool always_defringe)
{
  
  FitsImage *elixir = new FitsImage(Rim.ElixirName());
  if (!elixir->IsValid())
    {
      cerr << "elixir_to_toads : invalid file : "  << Rim.Name() << endl;
      return;
    }

  // copy and forget the input
  FitsImage image(Rim.FitsName(),*elixir,*elixir);
  delete elixir; elixir = NULL;


  // trim
  Frame Illu = TotalIlluRegion(image);
  image.Trim(Illu);
  

  // remove fringes
  string band = StringToUpper(image.KeyVal("TOADBAND"));
  if (band == "I" || band =="Z") {
    if(old_fringe_method) {
      if(! FringeUtils::IsDefringed(image) || always_defringe ) {
	cout << " removing fringes " << endl;
	FitsImage fringemap(Rim.FitsFringeName());
	if(!fringemap.IsValid()) {
	  cout << "ERROR in elixir_to_toads , fringemap " << Rim.FitsFringeName() << " is not valid" << endl;
	  return;
	}
	RemoveFringes(image, fringemap, true);
	image.AddOrModKey("FRINGED","SUB","Fringe pattern subtracted");
      }
    }else{
      if(! FringeUtils::IsDefringed(image) || always_defringe ) {
	if(! FringeUtils::IsDefringedWithNewMethod(image)) {
	  cout << " removing fringes with the new method" << endl;
	  float nsigma = 3;
	  if(band == "Z")
	    nsigma = 5;
	  if(FringeUtils::RemoveFringes(image,Rim.FitsFringeName(),0,nsigma,0,1)!=0) {
	    cout << "ERROR in elixir_to_toads at FringeUtils::RemoveFringes" << endl;
	  }
	  image.AddOrModKey("FRINGED","SUBNEW","Fringe pattern subtracted (new method)");
	}
      }
    }
  }

  // find saturation level
  float maxVal = image.MaxValue();   
  double saturation = maxVal, sat;
  double scale = 50.0;
  int loop = 0;
  while (saturation==maxVal && loop<2)
    {
      loop++;
      scale *= 10.0;
      Histo1d histo(int(maxVal/scale), 0, maxVal + 1 ); //HC
      Pixel *p = image.begin();
      Pixel *pend = image.end();
      for ( ; p < pend ; ++p) histo.Fill(*p, 1 ); // HC      
      double xMax;
      double maxContent = histo.MaxBin(xMax);
      double maxcoup = xMax;
      int count =1, bestMaxCount =0;
      for (int l=0; l<20; l++)
	{
	  histo.ZeroBins(maxVal*0.05*l, maxVal*0.05*(l+1));
	  sat = xMax;
	  maxContent = histo.MaxBin(xMax);
	  if (sat == xMax)
	    {
	      count++;
	      if (count >= bestMaxCount && sat != maxcoup)
		{ 
		  bestMaxCount = count;
		  saturation = sat;
		}
	    }
	  else count = 1;
	}
    }
  image.AddOrModKey("SATURLEV",saturation," Current saturation level");

  if (!bool(image.KeyVal("BACK_SUB")))
    {
      // compute stats
      Pixel sky,sig;
      image.SkyLevel(&sky, &sig);
      double rdnoise = image.KeyVal("TOADRDON");
      double gain = image.KeyVal("TOADGAIN");
      double sigth = sqrt(sky*gain+rdnoise*rdnoise)/gain;
      image.AddOrModKey("SKYLEV",sky," Original sky level");
      image.AddOrModKey("SKYSIGEX",sig," Original sky sigma");
      image.AddOrModKey("SKYSIGTH",sigth," app. Theoretical sigma (sqrt(sky*gain+rdnoise^2)/gain");
  
  // fix the masked 0 pixels with sky (assume non subtracted sky)
      Pixel *p = image.begin();
      Pixel *pend = image.end();
      for ( ; p < pend ; ++p) if (fabs(*p) < 0.1)  *p = sky;
    }
  else
    {
      cout << " this image does not look like a new born Elixir image : no sky statistics computed" << endl;
    }
}

void DumpHelp(const char *programname) {
  cout << programname << "  <DbImage(s)> " << endl 
       << " option:  --old_fringe_method (use old method to remove fringes)" << endl
       << "          --always_defringe (try to defringe image in any case (if the fringe pattern is a new one)" << endl
       << endl
       << "   will:\n" 
       << "     - trim the calibrated.fits \n"
       << "     - defringe if not yet defringed and i or z band \n"
       << "     - *NOT* multiply by the gain \n"
       << "     - compute the saturation and sky levels \n"
       << "     - fix the masked images if 0 \n"
       << endl;
  exit(-1);
}

int main(int argc,char **argv)
{
  if(argc <=1)
    DumpHelp(argv[0]);
  
  bool old_fringe_method = false;
  bool always_defringe = false;
  
  vector<string> filenames;
  for (int i=1; i< argc; ++i) {
    if (!strcmp(argv[i],"--old_fringe_method")) {
      old_fringe_method = true;
      continue;
    }
    if (!strcmp(argv[i],"--always_defringe")) {
      always_defringe = true;
      continue;
    }
    string filename = argv[i];
    filenames.push_back(filename);
  }
  
  int nimages = filenames.size();
  if(nimages<=0) {
    cout << "missing input image(s)" << endl;
    DumpHelp(argv[0]);
  }
  for(int im=0;im<nimages;im++) {
    ReducedImage im(filenames[im]);
    elixir_process(im,old_fringe_method,always_defringe);
  }

  return 1;
}

