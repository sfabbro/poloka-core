#include <string>

#include "fileutils.h"
#include "fitsimage.h"
#include "fitstoad.h"
#include "frame.h"
#include "reducedimage.h"
#include "histo1d.h"


static double FindSatur(const Image &image)
{
  
  // find saturation level
  // Why was this code copied from the library??
  // this is a shame!
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
  return saturation;
}


static void usage(const char * prog)
{
  cout << prog << " [-s (separate amps)] fitsimage ..." << endl;
  exit(1);
}

int main(int argc,char **argv)
{
  bool sep_amps = false;
  vector<string> filenames;
  for (int i=1; i< argc; ++i) {
    if (strcmp(argv[i],"-s") == 0)
      {
	sep_amps = true; continue;
      }
    string filename = argv[i];
    filenames.push_back(filename);
  }
  

  bool ok = true;
  size_t nimages = filenames.size();
  if(nimages<=0) {
    cout << "missing input image(s)" << endl;
    usage(argv[0]);
  }

  for(size_t i=0; i<nimages; i++) {
    FitsImage image(filenames[i]);
    int namp = image.KeyVal("TOADNAMP");
    cout << filenames[i] << ' ';
    if (sep_amps)
      for (int iamp=1; iamp<=namp; ++iamp)
	{
	  Frame frame=IlluRegion(image,iamp);
	  frame.xMax -= 1;
	  frame.yMax -= 1;
	  Image subimage=image.Subimage(frame);
	  cout << iamp << ' ' << FindSatur(subimage) << ' ';
	}
    else
      cout << FindSatur(image);
    cout << endl;
  }

  return EXIT_SUCCESS;
}

