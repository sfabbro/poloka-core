#include <iostream>
#include "fitsimage.h"

static void usage(char *progName)
{
  cerr << progName << " <LBL FITS files>" << endl;
}

int main(int nargs, char **args)
{
  if (nargs <=1) {usage(args[0]); exit(-1);}
  for (int i=1; i<nargs; ++i)
    {
      string name = args[i];
      FitsHeader imFits(name,RW);
      double satur = imFits.KeyVal("WELLDEPT");
      // LBL mulitplied every image by their gain but do not keep track 
      // of it anywhere in FITS files
      double sky = imFits.KeyVal("SKY");
#if 0 
      // old procedure trying to check if the WELLDEPT keyword was properly set
      // in LBL images. Useless now that most LBL images have been re-flatfielded properly.
      double gain = imFits.KeyVal("TOADGAIN");
      Pixel mean,sig;
      imFits.SkyLevel(&mean, &sig);//this is long, but nothing better as a check
      // check for sky subtracted : more than theoretical value measured on the frame
      if (fabs(mean-sig*sig)/mean > 5)
	{
	  imFits.AddOrModKey("BACK_SUB",true,"Sky was subtracted from LBL");
	}
      else sky = mean;
      if (satur < gain*40000) //hardcoded
	{
	  cout << " Changing saturation from " << satur << " to " << satur+sky << endl;
	  satur += sky;
	}
#endif
      imFits.AddOrModKey("SKYLEV",sky," Original sky level");
      imFits.AddOrModKey("SATURLEV",satur," Saturation (WELLDEPT (+sky in certain cases)");
      imFits.AddOrModKey("TOADGAIN",1,"Hopefully image is in photoelectrons");
    }
  return EXIT_SUCCESS;
}
