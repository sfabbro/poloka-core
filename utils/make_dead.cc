#include <iostream>
#include <vector>
#include "fitsimage.h"

void MakeDead(const Image &Flat, Image &Dead, 
	      const Pixel minv, const Pixel maxv)
{
  Pixel *pflat = Flat.begin();
  Pixel *pdead = Dead.begin();
  for (long i=Flat.Nx()*Flat.Ny(); i; --i, ++pflat, ++pdead)
    if ((*pflat < minv) || (*pflat > maxv)) *pdead = 1;
  
}

void NormalizeImage(Image &Im)
{
  Pixel mean,sigma;
  Im.SkyLevel(&mean, &sigma);
  cout << " norm = " << mean;
  Im /= mean;
}


static void usage(const char *progName)
{
  cerr << progName << " -l <value> -h <value> <FITS flat files>\n" 
       << "    -l : lowest value accepted by normalized flat  [0.6]\n"
       << "    -h : highest value accepted by normalized flat [1.4]\n";
}

int main(int nargs, char **args)
{
  if (nargs <=1) {usage(args[0]); exit(-1);}
  vector<string> toDo;
  Pixel minv=0.6,maxv=1.4;
  for (int i=1; i<nargs; ++i) 
    {
      char *arg = args[i];
      if (arg[0] != '-') 
	{
	  toDo.push_back(arg);
	  continue;
	}
      switch (arg[1])
	{
	case 'l': ++i; minv = atof(args[i]); break;
	case 'h': ++i; maxv = atof(args[i]); break;
	default : usage(args[0]); exit(1);
	}
    }

  int nim = toDo.size();

  cout << " " << nim << " images to process \n";
  for (int i=0; i<nim; ++i)
    {
      string name = toDo[i];
      FitsImage flat(name);
      cout << " Processing " << name;
      string deadname = "dead_"+name;
      FitsImage dead(deadname, (FitsHeader) flat);
      NormalizeImage(flat);
      MakeDead(flat, dead, minv, maxv);
      dead.ModKey("BITPIX",8);
      cout << endl;
    }
  return EXIT_SUCCESS;
}

