#include <iostream>
#include <string>
#include <dimage.h>

using namespace std;

int main (int argc, char ** argv ) {
  if(argc<2)
    return -1;
  
  for(int i=1;i<argc;i++) {
    string vignetfits = argv[i];
    Kernel toto;
    toto.readFits(vignetfits);
    int hx = toto.HSizeX()-4; // rm borders
    int hy = toto.HSizeY()-4;
    Kernel newtoto(hx,hy);
    for(int j=-hy;j<=hy;j++)
      for(int i=-hx;i<=hx;i++)
	if(toto(i,j)>0.01) // rm negative values
	  newtoto(i,j)=toto(i,j);
	else
	  newtoto(i,j)=0;
    newtoto.writeFits(vignetfits+"_4gif");
  }
  
  return 1;

}
