#include "fitsimage.h"
#include "dbimage.h"
#include "superflat.h"


static void usage(const char *progName)
{
  cout << endl 
       << progName << " removes fringes of a FITS image" << endl
       << "   usage : defringe <image in> -f <fringemap> [options]"<< endl
       << "       or  defringe <dbimage> [options]"<< endl
       << "    [options]:" << endl
       << "        -o <file> : output defringed image (default is overwrite)" << endl
       << "        -n        : normalize image to minimize rms sky" << endl
       << endl;
  exit(-1);
}



int main(int argc, char**argv)
{

  if (argc < 2) usage(argv[0]);

  string in_name = "", out_name = "", fringemap_name="";
  bool normflag = false;
  for (int i=1; i< argc; i++)
    {
      char *arg = argv[i];
      if (arg[0] != '-') {in_name=arg;}
      else switch(arg[1])
	{
	case 'o': i++; out_name=argv[i]; break;
	case 'f': i++; fringemap_name=argv[i]; break;
	case 'n': normflag = true;break;
	default : usage(argv[0]);
	}
    }

  FitsFileMode fmode = RW;
  if (out_name != "") fmode = RO;
  
  if (fringemap_name == "") 
    {
      DbImage dbim(in_name);
      in_name = dbim.FitsImageName(Calibrated);
      fringemap_name = dbim.FitsFringeName();
    }

  FitsImage to_fringe(in_name, fmode);
  FitsImage fringe_map(fringemap_name);
  Pixel sky,sigma;
  to_fringe.SkyLevel(&sky, &sigma);
  cout << " Before removing fringes : Sky = " << sky << "  Sigma = " << sigma << endl;

  RemoveFringes(to_fringe, fringe_map, normflag);

  to_fringe.SkyLevel(&sky, &sigma);
  cout << " After removing fringes : Sky = " << sky << "  Sigma = " << sigma << endl;

  if (fmode == RO)
    {
      FitsImage outFits(out_name, to_fringe, to_fringe);
      outFits.AddOrModKey("SKYLEV",sky,"Sky Level in photoelectrons");
      outFits.AddOrModKey("SKYSIGEX",sigma,"Sigma of Sky Level obtained from Image");
      outFits.AddOrModKey("FRINGED","SUB","Fringe pattern subtracted");
    }
  else 
    {
      to_fringe.AddOrModKey("SKYLEV",sky,"Sky Level in photoelectrons");
      to_fringe.AddOrModKey("SKYSIGEX",sigma,"Sigma of Sky Level obtained from Image");
      to_fringe.AddOrModKey("FRINGED","SUB","Fringe pattern subtracted");
    }

  return 1;
}
