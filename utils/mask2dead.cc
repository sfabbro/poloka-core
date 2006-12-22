#include <iostream>
#include <vector>

#include "fitsimage.h"
#include "fitstoad.h"
#include "frame.h"
#include "polokaexception.h"


using namespace std;

static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " deadinput.{fits,fz} ( deadoutput.{fits,fz} )  " << endl;
  cerr << " converts Megacam masks to poloka stuff: trims (if needed) inverts, and converts to BITPIX=8 (if output file != input file)" << endl;
  exit(-1);
}

int main(int nargs, char **args)
{
  // if nothing is given
  if (nargs < 2){usage(args[0]);}
  string inputname = args[1];
  string outputname = "";
  FitsFileMode  iomode = RW;
  if(nargs>2) {
    iomode = RO;
    outputname = args[2];
  }
  

  
  try{ 



  FitsImage deadinput(inputname,iomode);
  if ( ! deadinput.IsValid()) {
    cerr << "error reading " << args[1] << endl;
    usage(args[1]);
  }
  


  FitsImage *image=0;
  if(iomode==RW)
    image = &deadinput;
  else
    {
      image = new FitsImage(outputname,(const FitsHeader &)deadinput);
      image->ModKey("BITPIX",8);
      (Image &)(*image) = (const Image &) deadinput;
    }

  Frame illu = TotalIlluRegion(*image);
  image->Trim(illu);

  cout << image->Nx() << " " << image->Ny() << endl;
  
  for(int j=0;j<image->Ny();j++) {
    for(int i=0;i<image->Nx();i++) {
      if ((*image)(i,j)>0.5)
	(*image)(i,j)=0;
      else
	(*image)(i,j)=1;
    }
  }
  if (iomode != RW) delete image; // mandatory : it actually writes
  return EXIT_SUCCESS;



  }catch(PolokaException p) {
    p.PrintMessage(cout);
  }
  return EXIT_FAILURE;
}
