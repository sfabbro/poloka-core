#include <iostream>
#include <iomanip>
#include <string>
#include "fitsimage.h"
//#include "../dao_stuff/cdaophot.h"

int main(int argc,char **args)
{
  if (argc <=1)
    {
      cout << " image_stats <fitsfile(s)> " << endl;
      exit(-1);
    }
  size_t nchar = 10;
  for (int i=1; i< argc; ++i) if (strlen(args[i]) > nchar) nchar = strlen(args[i]);

  cout << endl << setiosflags(ios::left)
       << setw(nchar) << "IMAGE" 
       << setiosflags(ios::right)
       << setw(9) << "MEAN"
       << setw(7) << "SIG" 
       << setw(7) << "STH" 
       << setw(7) << "STHG1" 
       << setw(7) << "STHG" 
       << resetiosflags(ios::right)
       << endl << endl;

  for (int i=1; i< argc; ++i)
    {
      FitsImage image(args[i]);
      cout << setiosflags(ios::left) 
	   << setw(nchar) << args[i];

      if (!image.IsValid())
	{
	  cout << " invalid file "  << endl;
	  continue;
	}
      Pixel mean,sigma;
      image.SkyLevel(&mean, &sigma);

      float gain = image.KeyVal("TOADGAIN");
      float oldgain = image.KeyVal("GAIN");
      if (oldgain==0) oldgain = image.KeyVal("OLDGAIN");

      float rdnoise = image.KeyVal("TOADRDON");

      float sigth = sqrt(mean*gain+rdnoise*rdnoise)/gain;
      float sigth1 = sqrt(mean + rdnoise*rdnoise);
      float sigthold = sqrt(mean*oldgain+rdnoise*rdnoise)/oldgain;

      cout << setiosflags(ios::right) << setiosflags(ios::fixed)
	   << setw(9) << setprecision(2) << mean 
	   << setw(7) << setprecision(2) << sigma 
	   << setw(7) << setprecision(2) << sigth 
	   << setw(7) << setprecision(2) << sigth1 
	   << setw(7) << setprecision(2) << sigthold
	   << resetiosflags(ios::right)
	   << endl;
    }
  
  return 1; 
}
