#include <reducedimage.h>
#include <daophotpsf.h> 
#include <dimage.h> 




int main(int argc, char** argv) {

  if(argc<2)
    return -1;
  ReducedImage image(argv[1]);
  DaoPsf psf(image);
  
  Point Pt;
  int step = 40;
  int hy = step;
  int hx = step;
  
  DPixel ppdx,ppdy,psfval;
  DPixel integrale;

  // loop in all image area
  if(true) {
  for(int jc = step+5 ; jc < image.YSize()-step-5; jc += 100)
    for(int ic = step+5 ; ic < image.XSize()-step-5; ic += 100) {
      Pt.x = float(ic);
      Pt.y = float(jc);
      integrale = 0;
      for(int j=-hy;j<=hy;++j) {
	for(int i=-hx;i<=hx;++i) {
	  psfval = psf.Value(i+ic,j+jc, Pt, ppdx, ppdy);
	  integrale += psfval;
	  if( (!(psfval>0)) && (!(psfval<=0)) ) { // i.e. nan
	    //if(!(psfval>0) ) {
	    cout << ic << "," << jc << " " 
		 << i << "," << j << " "
		 << psfval << endl;
	    abort();
	  }
	}
      }
      cout <<  ic << "," << jc << " integral= " << integrale << endl;
    }
  }else{
  
  
  Pt.x = 1933.75;
  Pt.y = 3728.15;
  int ic = int(Pt.x);
  int jc = int(Pt.y);
  hx = 29;
  hy = 29;
  integrale = 0;
  for(int j=-hy;j<=hy;++j) {
  for(int i=-hx;i<=hx;++i) {
    //int i = -15;
    //int j = -29;
      psfval = psf.Value(i+ic,j+jc, Pt, ppdx, ppdy);
      integrale += psfval;
      if( (!(psfval>=0)) && (!(psfval<0)) ) { // i.e. nan
	//if(!(psfval>0) ) {
	cout << ic << "," << jc << " " 
	     << i << "," << j << " "
	     << psfval << endl;
	    abort();
      }
	      }
       }
  cout <<  ic << "," << jc << " integral= " << integrale << endl;
  }
  
  return 0;
}
