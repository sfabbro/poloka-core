#include <iostream>
#include <string>
#include <algorithm>

#include <fitsimage.h>
#include <Magick++.h> 

// #define DEBUG

void usage (const char* arg) {
  cerr << "usage: " << arg << " <fitsimage> <weightimage> <gifimage> " << endl;
  cerr << "options : " << endl;
  cerr << "  -d     : show dead instead of discarding" << endl;
  cerr << "  -g #   : gamma correction" << endl;
  cerr << "  -n     : negate" << endl;
  cerr << "  -e     : equalize" << endl;
  cerr << "  -l     : log scale" << endl;
  cerr << "  -r     : geometric ratio (default = 5)" << endl;
  exit(EXIT_FAILURE);
}

using namespace std;



double filter(const double &x,bool aslog) {
  if(aslog) {
    if(x<1.) 
      return 0.;
    else
      return log(x);
  }
  return x;
}

int main(int argc,char **argv) {
  
  if(argc<4)
    usage(argv[0]);
  
  bool negate    = false;
  bool equalize  = false;
  float gamma    = 0.;
  bool aslog     = false;
  bool show_dead = false;
  int geom_ratio = 5;
  string fitsimage_filename   = argv[1];
  string fitsweight_filename = argv[2];
  string gifimage_filename    = argv[3];
  for(int i=4;i<argc;i++) {
    const char *arg = argv[i];
    if(! arg[0]=='-')
      usage(argv[0]);
    switch(arg[1]) {
    case 'g' :
      gamma = atof(argv[++i]);
      break;
    case 'l' :
      aslog = true;
      break;
    case 'n' :
      negate = true;
      break;
    case 'e' :
      equalize = true;
      break;
    case 'd' :
      show_dead = true;
      break;
    case 'r' :
      geom_ratio = atoi(argv[++i]);
      if(geom_ratio<1) {
	cerr << "warning, force geom_ratio to 1" << endl;
	geom_ratio = 1;
      }
      break;
      
    default :
      usage(argv[0]);
    }
  }
  
  if(argc>4)
    gamma    = atof(argv[4]);
  
  float pixmin,pixmax;
  
  FitsImage image(fitsimage_filename.c_str());
  FitsImage weight(fitsweight_filename.c_str());
  
  int nx = image.Nx();
  int ny = image.Ny();
  double val;
  pixmin = 1.e30;
  pixmax = -1.e30;
  Pixel* image_pix  = image.begin();
  Pixel* weight_pix = weight.begin();
  for(; image_pix!=image.end(); ++image_pix , ++weight_pix ) {
    if(*weight_pix) {
      val = filter(*image_pix,aslog);
      if (val < pixmin ) pixmin = val;
      else if (val > pixmax ) pixmax = val;
    }
  }
  
  
  int gif_nx = nx/geom_ratio;
  int gif_ny = ny/geom_ratio;
  
  Magick::Image  magick_image(Magick::Geometry(gif_nx,gif_ny) , "black" );
  magick_image.magick("gif"); 
  magick_image.type(Magick::TrueColorType);
  magick_image.modifyImage(); 
  Magick::PixelPacket*  colorpixels 
    = magick_image.getPixels(0,0,gif_nx,gif_ny);
  
  for(int j=0;j<gif_ny ; j++) {
    for(int i=0;i<gif_nx ; i++ , colorpixels++) {
      double sum_w = 0;
      double sum_wf = 0;
      bool dead = false;
      for(int l=j*geom_ratio;l<(j+1)*geom_ratio;l++) {
	for(int k=i*geom_ratio;k<(i+1)*geom_ratio;k++) {
	  const double &w = weight(k,l);
	  dead = (w==0 && show_dead);
	  if(dead) break;
	  sum_w += w;
	  sum_wf += w*filter(image(k,l),aslog);
	}
	if(dead) break;
      }
      if(dead)
	val = 0;
      else
	val = (sum_wf/sum_w - pixmin)/(pixmax-pixmin); 
      *colorpixels = Magick::ColorRGB(val,val,val);
    }
  }
  
  try {
    //view.sync();
    magick_image.syncPixels();
  }catch( Magick::Exception &error_ ) {
    cerr << "Caught exception when synchronizing images : " << error_.what() << endl;
    return 1;
  }
  
  if(gamma>0.) magick_image.gamma(gamma);
  if(negate)   magick_image.negate();
  if(equalize) magick_image.equalize();
  
  magick_image.write(gifimage_filename.c_str());
  
  return 0;
}
