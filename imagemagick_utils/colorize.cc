#include <iostream>
#include <string>
#include <algorithm>

#include <fitsimage.h>
#include <Magick++.h> 

// #define DEBUG

void usage (const char* arg) {
  cerr << "usage: " << arg << " -i <fitsimage1> <color1> <scale1> ( -i <fitsimage2> <color2> <scale2> ....)  -o <gifimage>" << endl;
  cerr << "options :  -m <minvalue>" << endl;
  cerr << "           -M <maxvalue>" << endl;
  cerr << "           -n : negate" << endl;
  exit(EXIT_FAILURE);
}

using namespace std;

void addfitsinmagick(const Image& fitsimage,Magick::Image &magickimage, const Magick::ColorRGB &color, float pixmin, float pixmax, float scale=1) {
  // color 
  double blue_toadd  = color.blue();
  double green_toadd = color.green();
  double red_toadd  = color.red();
  
  // deal with geometry
  unsigned int fits_nx = fitsimage.Nx();
  unsigned int fits_ny = fitsimage.Ny();
  Magick::Geometry geom;
  try {
    geom = magickimage.size();
  }catch( Magick::Exception &error_ ) {
    cerr << "Caught exception when reading geometry: " << error_.what() << endl;
    return;
  }
  unsigned int magick_nx = geom.width();
  unsigned int magick_ny = geom.height();  
  if(fits_nx > magick_nx || fits_ny>magick_ny ) {
    cerr << "in addfitsinmagick, magick image is smaller than fits image, cannot handle this" << endl;
    return;
  }
  // if fitsimage is smaller, add it at the center (tuned for kernel display)
  
  //Magick::Pixels view(magickimage);
  Magick::PixelPacket*  colorpixels 
    = magickimage.getPixels((magick_nx-fits_nx)/2,(magick_ny-fits_ny)/2,fits_nx,fits_ny); // to modify
  
  float pixval,scaledval;
  for(unsigned int j=0;j<fits_ny;j++) {
    for(unsigned int i=0;i<fits_nx;i++) {
      pixval = fitsimage(i,j)*scale;
      scaledval = (pixval - pixmin)/(pixmax - pixmin);
      if(scaledval<0) scaledval=0;
      if(scaledval>1) scaledval=1;
      
      Magick::ColorRGB pixcolor = Magick::Color(*colorpixels);
#ifdef DEBUG
      cout << " init " 
	   << pixcolor.blue() << " "
	   << pixcolor.green() << " "
	   << pixcolor.red() << " "
	   << endl;
      
      cout << i << " " 
	   << j << " " 
	   << scaledval  << " "
	   << blue_toadd << " "
	   << green_toadd << " "
	   << red_toadd << " "
	   << endl;
#endif
     
      pixcolor.blue(min(pixcolor.blue()+blue_toadd*scaledval,1.0));
      pixcolor.green(min(pixcolor.green()+green_toadd*scaledval,1.0));
      pixcolor.red(min(pixcolor.red()+red_toadd*scaledval,1.0));
      *colorpixels++ = pixcolor;
      
#ifdef DEBUG
      cout << " end " 
	   << pixcolor.blue() << " "
	   << pixcolor.green() << " "
	   << pixcolor.red() << " "
	   << endl;
#endif
    }
  }
  
  try {
    //view.sync();
    magickimage.syncPixels();
  }catch( Magick::Exception &error_ ) {
    cerr << "Caught exception when synchronizing images : " << error_.what() << endl;
    return;
  }

#ifdef DEBUG
  cout << "===================================" << endl;
#endif
}

void getmaxsizeandscale(const vector<string> &fitsfilenames,int &nx,int &ny, float &pixmin, float &pixmax) {
  nx=ny=0;
  pixmin=100000;
  pixmax=-100000;
  for(unsigned int i=0;i<fitsfilenames.size();i++) {
    FitsImage img(fitsfilenames[i].c_str());
    int im_nx = img.Nx();
    int im_ny = img.Ny();
    float im_pixmin,im_pixmax;
    img.MinMaxValue(&im_pixmin,&im_pixmax);
    if(nx<im_nx)
      nx=im_nx;
    if(ny<im_ny)
      ny=im_ny;
    if(im_pixmin<pixmin)
      pixmin=im_pixmin;
    if(im_pixmax>pixmax)
      pixmax=im_pixmax;
  }
}

Magick::ColorRGB colorfromname(const string &colorname) {
  char colorid = colorname[0];
  Magick::ColorRGB color;
    switch(colorid) {
    case 'w':
    case 'W':
      color.blue(1);
      color.green(1);
      color.red(1);
      break;
    case 'g':
    case 'G':

      color.blue(1);
      color.green(0.1);
      color.red(0);
      break;
    case 'r':
    case 'R':
      color.blue(0.1);
      color.green(1); 
      color.red(0.1);
      break;
    case 'i':
    case 'I':
      color.blue(0);
      color.green(0.1);
      color.red(1);
      break;
    case 'z':
    case 'Z':
      color.blue(0);
      color.green(0);
      color.red(1);
      break;
    default:
      cerr << "unknown color " << colorid << endl;
      usage("ERROR");
    }
    return color;
}

int main(int argc,char **argv) {
  
  if(argc<4)
    usage(argv[0]);
  
  int i=0;
  string outputfilename="color.gif";
  vector<string> fitsfilenames;
  vector<Magick::ColorRGB> colors;
  vector<float> scales;
  float pixmin=10000;
  float pixmax=-10000;
  bool negate = false;
  for (int i=1; i<argc; ++i) {
    char *arg = argv[i];
    if (arg[0] != '-') {
      cerr << "unexpected argument " << arg << endl;
      usage(argv[0]);
    }
    switch (arg[1]) {
    case 'n':
      negate = true;
      break;
    case 'i': 
      fitsfilenames.push_back(argv[++i]);
      colors.push_back(colorfromname(argv[++i]));
      scales.push_back(atof(argv[++i]));
      break;
    case 'M':
      pixmax = atof(argv[++i]);
      break;
    case 'm':
      pixmin = atof(argv[++i]);
      break;
    case 'o':
      outputfilename = argv[++i];
      break;
    default:
      cerr << "unknown option " << arg << endl;
      usage(argv[0]);
      break;     
    }
  }
  
  
  unsigned int nimages=fitsfilenames.size();
  
  int nx,ny;
  float lpixmin,lpixmax;
  getmaxsizeandscale(fitsfilenames,nx,ny,lpixmin,lpixmax);
  
  if(pixmin>9000)
    pixmin=0; // force min pix to 0
  if(pixmax<-9000)
    pixmax=lpixmax;
  
  cout << "MinMax " << pixmin << " " << pixmax << endl;
  
  Magick::Image  magic_image(Magick::Geometry(nx,ny) , "black" );
  magic_image.magick("gif"); 
  magic_image.type(Magick::TrueColorType);
  magic_image.modifyImage(); 
  
  for(unsigned int i=0;i<nimages;i++) {
    FitsImage fits_image(fitsfilenames[i].c_str());
    addfitsinmagick(fits_image,magic_image,colors[i],pixmin,pixmax,scales[i]);
  }
 
  /*
  int scale = 4;
  Magick::Geometry newgeom(nx*scale,ny*scale);
  magic_image.scale(newgeom);
  */
  if(negate)
    magic_image.negate();
  magic_image.gamma(3.); 
  //magic_image.normalize();

  magic_image.write(outputfilename.c_str());
  
  return 0;
}
