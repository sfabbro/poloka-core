#include <iostream>

#include <fitsimage.h>

using namespace std;

int usage(char* pg) {
  cout << pg << ": performs very basic aperture photometry" << endl; 
  cout << "usage : " << pg << " <image> x y radius (pixels) " << endl;
  return 0;
}

int main(int argc, char** argv) {
  if(argc !=5) {
    return usage(argv[0]);
  }
  
  FitsImage image(argv[1]);
  if(!image.IsValid()) {
    cerr << "image " << argv[1] << "is not valid" << endl;
    return usage(argv[0]);
  }
  float x = atof(argv[2]);
  float y = atof(argv[3]);
  float radius = atof(argv[4]);
  
  double sum=0;
  for(int i=x-radius;i<x+radius;i++) {
    if(i<0) continue;
    if(i>=image.Nx())
      break;
    for(int j=y-radius;j<y+radius;j++) {
      if(j<0) continue;
      if(j>=image.Ny())
	break;
      sum += image(i,j);
    }
  }
  
  cout << sum << endl;
      


  return 0;
}
