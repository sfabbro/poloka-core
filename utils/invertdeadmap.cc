#include<fitsimage.h>
#include "polokaexception.h"

int main (int argc , char **argv) {

  if(argc<2) {
    cout << argv[0] << " <fitsimage> (dead map to invert)" << endl;
    exit(-12);
  }

  bool ok = true;
  for(int i=1;i<argc;++i) {
    
    try{

      FitsImage tutu(argv[i],RW);
      
      Pixel mean,rms;
      tutu.SkyLevel(&mean,&rms);
      if(mean>0.5) {
	cout << "inverting dead map " << argv[i] << endl;
	tutu *= -1.;
	tutu += 1;
      }else{
	cout << "dead map " << argv[i] << " is ok" << endl;
      }
      
    }catch(PolokaException p) {
      p.PrintMessage(cout);
      ok=false;
    }
    
  }
  
  if(ok)
    return EXIT_SUCCESS;
  return EXIT_FAILURE;
  
}


