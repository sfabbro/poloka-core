#include <iostream>
#include <list>
#include <poloka/fitsimage.h>

static void MakeDead(const Image &Flat, Image &Dead, 
	      const Pixel& minv, const Pixel& maxv)
{
  Pixel *pflat = Flat.begin();
  Pixel *pdead = Dead.begin();
  for (long i=Flat.Nx()*Flat.Ny(); i; --i, ++pflat, ++pdead)
    if ((*pflat < minv) || (*pflat > maxv)) *pdead = 1;
  
}

static void NormalizeImage(Image &Im)
{
  Pixel mean,sigma;
  Im.SkyLevel(&mean, &sigma);
  Im /= mean;
}


static void usage(const char *progname)
{
  cerr << "Usage: " << progname << "[OPTIONS]... FITS...\n"
       << "Create a dead pixel map from a FITS flat file\n\n"
       << "    -l FLOAT : lowest value accepted by normalized flat (default: 0.6)\n"
       << "    -h FLOAT : highest value accepted by normalized flat (default: 1.4)\n";
  exit(EXIT_FAILURE);
}

int main(int nargs, char **args)
{
  if (nargs <=1) usage(args[0]);

  list<string> imList;
  Pixel minv=0.6,maxv=1.4;
  for (int i=1; i<nargs; ++i) 
    {
      char *arg = args[i];
      if (arg[0] != '-') 
	{
	  imList.push_back(arg);
	  continue;
	}
      switch (arg[1])
	{
	case 'l': ++i; minv = atof(args[i]); break;
	case 'h': ++i; maxv = atof(args[i]); break;
	default : usage(args[0]);
	}
    }

  for (list<string>::iterator it = imList.begin(); it != imList.end(); ++it)
    {
      string name = *it;
      FitsImage flat(name);
      cout << " Produce dead pixel map for " << name << endl;
      string deadname = "dead_"+name;
      FitsImage dead(deadname, (const FitsHeader &) flat);
      NormalizeImage(flat);
      MakeDead(flat, dead, minv, maxv);
      dead.ModKey("BITPIX",8);
    }
  return EXIT_SUCCESS;
}
