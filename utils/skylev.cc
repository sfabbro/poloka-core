#include <iostream>
#include <string>

#include "fitsimage.h"
#include "frame.h"
#include "fitstoad.h"
#include "alltelinst.h"
#include "polokaexception.h"


static void usage(const char *prog)
{
  cout << prog << " <fitsfile(s)> " << endl
       << " [-s (to separate levels for diff. amplifiers)] " 
       << endl;
    exit(EXIT_FAILURE);
}  

int main(int argc,char **args)
{
  bool splitAmps = false;
  bool ok = true;
  if (argc <=1)
    {
      usage(args[0]);
    }
  for (int i=1; i< argc; ++i)
    {
      const char *arg = args[i];
      if (arg[0] == '-')
	{
	  switch (arg[1])
	    {
	    case 's' : splitAmps = true; break;
	    default : usage(args[0]);
	    }
	  continue;
	}
      try
	{
	  FitsImage image(arg);
	  int Namp;
	  if (splitAmps && image.HasKey("TOADNAMP")) Namp = image.KeyVal("TOADNAMP");
	  else Namp = 1;
      
	  Pixel mean,sigma;
	  if (Namp==1)
	    {
	      FitsHeader  &head = image;
	      Pixel meanOS = 0, sigmaOS =0;
	      if (!head.HasKey("WRITEDAT"))
		{
		  Frame overscan= OverscanRegion(head,Namp);
		  if (overscan.Area() != 0)
		    {
		      cout << " raw image: subtract bias from the mean" << endl;
		      meanOS = image.MedianInFrame(overscan, sigmaOS);
		      cout << "Median of overscan region " << meanOS << endl;
		    }
		}
	      image.SkyLevel((head),& mean, &sigma);
	      cout << args[i] << ' ' << "  m : " << mean - meanOS<< " s : " << sigma << endl;
	    }
	  else // several amps
	    {
	      cout << args[i] << " : " << endl;
	      Frame TotalIllu = TotalIlluRegion(image);
	      for (int k=1; k<= Namp; ++k)
		{
		  Frame region;
		  if (image.Nx() == TotalIllu.Nx() && image.Ny() == TotalIllu.Ny()) region = AmpRegion(image,k);
		  else region = IlluRegion(image,k);
		  image.SkyLevel(region,&mean,&sigma);
		  cout << " ampli("<< k  <<")  m : " << mean << " s : " << sigma << endl;
		}
	    }
	}
      catch(PolokaException p)
	{
	  p.PrintMessage(cout);
	  ok = false;	      
	} 
    }// end of loop on args
  return ((ok == true) ? EXIT_SUCCESS : EXIT_FAILURE);

}
