#include <iostream>
#include <vector>

#include "fitsimage.h"


using namespace std;

static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " deadinput.{fits,fz} ( deadoutput.{fits,fz} )  " << endl;
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
  
  
  FitsImage deadinput(inputname,iomode);
  if ( ! deadinput.IsValid()) {
    cerr << "error reading " << args[1] << endl;
    usage(args[1]);
  }
  
  Image *image=0;
  if(iomode==RW)
    image = (Image*)(&deadinput);
  else
    image = new Image(deadinput.Nx(),deadinput.Ny());
  
  cout << image->Nx() << " " << image->Ny() << endl;
  for(int j=0;j<image->Ny();j++) {
    for(int i=0;i<image->Nx();i++) {
      if(deadinput(i,j)>0.5)
	(*image)(i,j)=0;
      else
	(*image)(i,j)=1;
    }
  }
  if(iomode==RO) {
    FitsImage toto(outputname,deadinput,*image);
    toto.PreserveZeros();
  }else{
    deadinput.PreserveZeros(); 
  }
  return 0;
}
