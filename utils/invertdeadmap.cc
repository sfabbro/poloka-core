#include<fitsimage.h>

int main (int argc , char **argv) {

  if(argc<2) {
    cout << argv[0] << " <fitsimage> (dead map to invert)" << endl;
    exit(-12);
  }

  for(int i=1;i<argc;++i) {
    
    FitsImage tutu(argv[i],RW);
    if(!tutu.IsValid()) {
      cout << argv[1] << "is pas cool" << endl;
      exit(-40);
    }
    
    Pixel mean,rms;
    tutu.SkyLevel(&mean,&rms);
    if(mean>0.5) {
      cout << "inverting dead map " << argv[i] << endl;
      tutu *= -1.;
      tutu += 1;
    }else{
      cout << "dead map " << argv[i] << " is ok" << endl;
    }
  }
  return 0;
}


