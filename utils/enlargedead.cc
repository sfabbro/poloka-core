#include <iostream>
#include <vector>

#include "fitsimage.h"


using namespace std;

static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " deadinput.{fits,fz} deadoutput.{fits,fz}  " << endl;
  exit(-1);
}

int main(int nargs, char **args)
{
  // if nothing is given
  if (nargs < 3){usage(args[0]);}
  
  
  FitsImage deadinput(args[1]);
  if ( ! deadinput.IsValid()) {
    cerr << "error reading " << args[1] << endl;
    usage(args[1]);
  }
  
  Image deadoutput(deadinput.Nx(),deadinput.Ny());
  int i2,j2;
  for(int j=0;j<deadoutput.Ny();j++) {
    j2 = (j/2)*2;
    for(int i=0;i<deadoutput.Nx();i++) {
      if(deadinput(i,j)>0.5) {
	i2 = (i/2)*2;
	deadoutput(i2,j2)=1;
	deadoutput(i2+1,j2)=1;
	deadoutput(i2+1,j2+1)=1;
	deadoutput(i2,j2+1)=1;
      }
    }
  }
  FitsImage toto(args[2],deadinput,deadoutput);
  toto.PreserveZeros();

  return 0;
}
