#include <iostream>
#include <iomanip>
#include <string>

#include <poloka/fitsimage.h>

static void usage(const char* progname) {
  cerr << "Usage: " << progname << " FITS...\n"
	 << "Compute sky, r.m.s., expected r.m.s., min and max values of FITS image\n";
  exit(EXIT_FAILURE);
}

int main(int argc,char **args)
{
  if (argc<=1) usage(args[0]);

  size_t nchar = 10;
  for (int i=1; i< argc; ++i) if (strlen(args[i]) > nchar) nchar = strlen(args[i]);

  cout << endl << setiosflags(ios::left)
       << setw(nchar) << "IMAGE" 
       << setiosflags(ios::right)
       << setw(9) << "SKY"
       << setw(7) << "SIG" 
       << setw(7) << "STH" 
       << setw(9) << "MIN" 
       << setw(9) << "MAX" 
       << resetiosflags(ios::right)
       << endl << endl;

  for (int i=1; i< argc; ++i)
    {
      FitsImage image(args[i]);
      cout << setiosflags(ios::left) 
	   << setw(nchar) << args[i];

      if (!image.IsValid())
	{
	  cerr << args[0] << ": " << args[i] << ": invalid file\n";
	  continue;
	}
      Pixel mean,sigma,minv,maxv;
      image.SkyLevel(&mean, &sigma);
      image.MinMaxValue(&minv, &maxv);

      float gain = image.KeyVal("TOADGAIN");
      float rdnoise = image.KeyVal("TOADRDON");
      float sigth = 0;
      if (gain>0 && rdnoise>0)
	sigth = sqrt(mean*gain + rdnoise*rdnoise) / gain;

      cout << setiosflags(ios::right) << setiosflags(ios::fixed)
	   << setw(9) << setprecision(2) << mean << ' '
	   << setw(7) << setprecision(2) << sigma << ' '
	   << setw(7) << setprecision(2) << sigth << ' '
	   << setw(9) << setprecision(1) << minv << ' '
	   << setw(9) << setprecision(1) << maxv
	   << resetiosflags(ios::right)
	   << endl;
    }
  
  return EXIT_SUCCESS; 
}
