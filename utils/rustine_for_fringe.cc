#include <iostream>
#include <vector>
#include <fstream>

#include "fitsimage.h"
#include "histo1d.h"


class mypixel{
public:
  mypixel(int ix,int iy,Pixel ival) {
    x=ix;
    y=iy;
    flux=ival;
  }
  int x;
  int y;
  Pixel flux;
};


static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " fitsimage deadimage  " << endl;
  exit(-1);
}

// fill a vector of mypixels with values < min or > max
int findbadpixels(std::vector<mypixel> &badpixels, const Image &image,  Pixel min, Pixel max) {
  Pixel value;
  for(int y=0;y<image.Ny();++y)
  //for(int y=2;y<image.Ny()-2;++y) // warning, we don't take into account top and bottom rows (2 pixels width)
    for(int x=0;x<image.Nx();++x) {
      value = image(x,y);
      if(value<min || value > max) {
	mypixel badpix(x,y,value);
	badpixels.push_back(badpix);
      }
    }
  return badpixels.size();
}

// fill a new image 0 or 1 if( values < min or > max)
int findbadpixels(Image &badimage, const Image &image,  Pixel min, Pixel max) {
  if(badimage.Nx()!= image.Nx() || badimage.Ny()!= image.Ny()) {
    cerr << "images need to have the same size" << endl;
    return 0;
  }
    
  
  Pixel *pix = badimage.begin();
  Pixel *end = badimage.end();
  for(;pix!=end;++pix) *pix=0;
  
    
  Pixel value;
  int nbad=0;
  for(int y=0;y<image.Ny();++y)
    //for(int y=2;y<image.Ny()-2;++y)// warning, we don't take into account top and bottom rows (2 pixels width)
    for(int x=0;x<image.Nx();++x) {
      value = image(x,y);
      if(value<min || value > max) {
	badimage(x,y)=1;
	nbad++;
      }
    }
  return nbad;
}



void getnpix(const Image& image, Pixel mean, Pixel thres, int &nlt , int &ngt, Pixel &min, Pixel &max) {
  Pixel *pix = image.begin();
  Pixel *end = image.end();
  min = 1.e30;
  max = 1.e-30;
  nlt=0;
  ngt=0;
  for(; pix != end; ++pix) {
    if(*pix>mean+thres) ngt++;
    if(*pix<mean-thres) nlt++;
    if(*pix>max) max = *pix;
    if(*pix<min) min = *pix;
  }
}

void fillhisto(Histo1d &histo, const Image &image) {
  Pixel *pix = image.begin();
  Pixel *end = image.end();
  Pixel hmin = histo.Minx();
  Pixel hmax = histo.Minx()+histo.Nx()/histo.Scale();
  int nout = 0;
  for(; pix != end; ++pix) {
    histo.Fill(*pix);
    if(*pix<hmin || *pix>hmax)
      nout ++;
  }
  cout << "npix ut of histo = " << nout << endl;
}

// first and last bin positions with zero entries when starting from the middle
void findedges_with_non_zero(const Histo1d &histo, Pixel mean,Pixel &min, Pixel &max) {
  int middlebin = histo.BinAt(mean);
  cout << "middlebin nx" << middlebin << " " << histo.Nx() << endl;
  const float* harray = histo.array();
  for(int bin = middlebin; bin<histo.Nx();++bin) {
    if(harray[bin]==0) {
      cout << "max " << bin << " " << histo.BinCenter(bin) << endl;
      max = histo.BinCenter(bin);
      break;
    }
  }
  for(int bin = middlebin; bin>=0;--bin) {
    if(harray[bin]==0) {
      cout << "min " << bin << " " << histo.BinCenter(bin) << endl;
      min = histo.BinCenter(bin);
      break;
    }
  }
}

// first and last bin positions with zero entries when starting from the middle
void findedges_with_derivative(const Histo1d &histo,Pixel &max) {
  
  const float* harray = histo.array();
  
  // get the bin for which content = max/100 (starting from the end (right) of the histo
  double maxofhisto = 0;
  for(int bin = histo.Nx()-1; bin>=0;--bin) {
    if(harray[bin]>maxofhisto) maxofhisto = harray[bin];
  }
  cout << "maxofhisto = " << maxofhisto << endl;
  
  double value2reach = maxofhisto/100.;
  int startatbin = 0;
  for(int bin = histo.Nx()-1; bin>=0;--bin) {
    if(harray[bin]>value2reach) {
      startatbin = bin;
      break;
    }
  }
  
  // compute most negative derivative starting from staratbin
  for(int bin = startatbin; bin<histo.Nx()-1;++bin) {
    if(harray[bin+1]==0 || harray[bin+1]>harray[bin]) {
      max = histo.BinCenter(bin+2); // +2 is tuning
      break;
    }
  }
  cout << " startatbin = " << startatbin << endl;
  cout << " max = " << max << endl;
}

int interpolate_badpixels(Image &interpolatedimage,const Image &badimage,Image &deadimage, const Image &inputimage) {
  //assume all sizes are ok, too lazy to check now
  // interpolatedimage and deadimage will be modified


  int nx = badimage.Nx();
  int ny = badimage.Ny();

  //Pixel *interpolatedpix = interpolatedimage.begin();
  //Pixel *badpix = badimage.begin();
  
  int hxp,hxm;
  int hyp,hym;
  int dmin;
  int nsuperbad = 0;
  int n_vertically = 0;
  int n_horizontally = 0;
  bool direction_is_found;
  bool use_horizontally;
  bool use_vertically;
  for(int y=0;y<badimage.Ny();++y)
  //for(int y=2;y<ny-2;++y) //we don't take into account top and bottom of images which are at 0
    for(int x=0;x<nx;++x) {
      if(badimage(x,y)) {
	// init
	dmin = 10000;
	direction_is_found = false;
	use_horizontally = false;
	// start interpolation
	// horizontally
	
	hxp=x;hxm=x;
	while(hxp<nx && (badimage(hxp,y))) {hxp++;} // continue while pix is bad to the right
	if(hxp<nx) { // if not can't interpolate horizontally
	  while(hxm>=0  && (badimage(hxm,y))) {hxm--;} // continue while pix is bad to the left
	  if(hxm>=0){ // it's ok
	    use_horizontally = true;
	    dmin = hxp-hxm; // horzontal distance of interpolation
	    if(dmin<4) { // that's ok (this is to go faster with dead columns)
	      //cout << "horizontally is well enough" << endl;
	      direction_is_found = true;
	    }
	  }
	}
	if(!direction_is_found) {
	  // vertically
	  hyp=y;hym=y;
	  while(hyp<ny && (badimage(x,hyp))) {hyp++;} // continue while pix is bad to the right
	  if(hyp<ny) { // if not can't interpolate vertically
	    while(hym>=0  && (badimage(x,hym))) {hym--;} // continue while pix is bad to the left
	    if(hym>=0){ // it's ok
	      if(hyp-hym<dmin) {
		dmin = hyp-hym;
		use_horizontally = false;
		use_vertically = true;
	      }
	    }
	  }
	}
	if(dmin>10) {
	  if(deadimage(x,y)<0.1) {
	    badimage(x,y)=3; // it is not dead and not interpolated: pas cool
	    deadimage(x,y)=1; // we set this to dead
	    nsuperbad++;
	    cout << "could not interpolate and not dead " << x << " " << y << endl;
	  }else{
	    badimage(x,y)=2; // it is dead and not interpolated
	  }
	  interpolatedimage(x,y) = 0;
	}else{
	  // now do the interpolation
	  if(use_horizontally) {
	    n_horizontally ++;
	    interpolatedimage(x,y) = (inputimage(hxp,y)*(x-hxm)+inputimage(hxm,y)*(hxp-x))/(hxp-hxm);
	    badimage(x,y)=-1;
	    //cout << "hor = " << interpolatedimage(x,y) << " " << (hxp-hxm) << endl;
	  } else if(use_vertically) {
	    n_vertically ++;
	    interpolatedimage(x,y) = (inputimage(x,hyp)*(y-hym)+inputimage(x,hym)*(hyp-y))/(hyp-hym);
	    badimage(x,y)=-2;
	    //cout << "ver = " << interpolatedimage(x,y) << " " << (hyp-hym) << endl;
	  } else {
	    cout << "i should not be here" << endl;
	    interpolatedimage(x,y) = 0;
	    nsuperbad++;
	  }
	}
      }
    }
  cout << "nsuperbad = " << nsuperbad << endl;
  cout << "n_vertically = " << n_vertically  << endl;
  cout << "n_horizontally = " << n_horizontally  << endl;
  
  return nsuperbad;
}


int main(int nargs, char **args)
{
  if(nargs<3)
    usage(args[0]);
  
  Image deadimage; // will be modified 
  { // open, copy and close
    FitsImage deadfitsimage(args[2]);
    deadimage = deadfitsimage;
  }
  // set 0.00001 to 0 and 0.9999 to 1 (maybe pbs in the fitsimage (due to BSCALE)
  {
    Pixel* pix = deadimage.begin();
    Pixel* pixend = deadimage.end();
    for(;pix!=pixend;++pix)
      if(*pix>0.5)
	*pix=1;
      else
	*pix=0;
  }
  
  
  FitsImage image(args[1]);
  
  Pixel mean,sigma;
  float edge_min,edge_max;
  image.SkyLevel(&mean,&sigma);
  cout << "mean  = " << mean  << endl;
  cout << "sigma = " << sigma << endl;
  
  
  Pixel hmin = mean - 10*sigma;
  Pixel hmax = mean + 10*sigma;
  int nbins = int((hmax-hmin)/(sigma/10.));
  Histo1d histo(nbins,hmin,hmax);
  fillhisto(histo,image);
  findedges_with_non_zero(histo,mean,edge_min,edge_max);
  findedges_with_derivative(histo,edge_max);
  
  // get number of pixels higher than edge_max
  float nsuspicious = 0;
  const float* harray = histo.array();
  for(int i=histo.BinAt(edge_max);i<histo.Nx();i++) {
    nsuspicious += harray[i];
  }
  
  cout << "edge_min   = " << edge_min  << endl;
  cout << "edge_max   = " << edge_max << endl;
  cout << "naftermax   = " << nsuspicious << endl;
  // std::vector<mypixel> badpixels;
  // int nbad = findbadpixels(badpixels,image,edge_min,edge_max)
  //int nbad = findbadpixels_and_interpolate(interpolatedimage,image,edge_min,edge_max);
  Image badimage(image.Nx(),image.Ny());
  int nbad = findbadpixels(badimage,image,edge_min,edge_max);
  cout << "nbad  = " << nbad << endl;

  if(nbad>image.Ny()*100) {
    cout << "nbad is too large, let's exit" << endl;
    string name = args[1]; name +="_histo";
    ofstream st(name.c_str());
    st << histo << endl;
    st.close();
    exit(-1);
  }

  Image interpolatedimage = image;
  int nsuperbad = interpolate_badpixels(interpolatedimage,badimage,deadimage,image);
  
  //badimage += 10*deadimage;
  
  // if pixel>0, it was not interpolated but fortunatly it's dead, so ...
//   Pixel *pix = badimage.begin();
//   Pixel *endpix = badimage.end();
//   int notinterpolated_but_dead = 0;
//   for(;pix!=endpix;++pix)
//     if(*pix>10)
//       notinterpolated_but_dead ++;
//   cout << "nsuperbad = " << nsuperbad << endl;
//   cout << "notinterpolated_but_dead = " << notinterpolated_but_dead << endl;
//   int nsuperbad2 = nsuperbad- notinterpolated_but_dead;
//   cout << "nsuperbad2 = " << nsuperbad-notinterpolated_but_dead << endl;

  // save cleaned image
  {
    string name = args[1];
    if(nsuperbad<150) // this is quite arbitrary
      name+="_clean";
    else
      name+="_clean_but_bad";
    FitsImage clean_image(name,image,interpolatedimage);
    char comment[300];
    clean_image.AddCommentLine("this image was processed by rustine_for_fringe");
    sprintf(comment,"nbad=%d, nsuperbad=%d, naftermax=%d, edge_min=%f, edge_max=%f",nbad,nsuperbad,int(nsuspicious),(edge_min),(edge_max));
    clean_image.AddCommentLine(comment);
    // image written when destroyed
  }
  
  // save image of bad pixels for debugging
  {string name = args[1]; name +="_bad.gz"; FitsImage bad(name,image,badimage);}
  
  // save image of dead pixels, this is quite important after for defringing
  {string name = args[1]; name +="_dead.gz"; FitsImage dead(name,image,deadimage); dead.ModKey("BITPIX",8);}
  

  // and save histogram
  {
    string name = args[1]; name +="_histo";
    ofstream st(name.c_str());
    st << histo << endl;
    st.close();
  }
  return 0;
}
