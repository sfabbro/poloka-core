// -*- C++ -*-
// 
// 
#include <getopt.h>
#include <assert.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

#include "fitsimage.h"
#include "vutils.h"


using namespace std;


void usage()
{
  cerr << "usage: checkelixirflats [OPTIONS]" << endl;
  cerr << "OPTIONS:" << endl;
  cerr << " -r <rawflat>       specify the raw flat" << endl;
  cerr << " -s <scatterflat>   specify the scatter flat" << endl;
  cerr << " -p <photflat>      specify the photometric flat" << endl;
  cerr << " -o <imgname>       img=p/(r*s)" << endl;
  cerr << " -n <ntname>        values" << endl;
  cerr << " -h                 print this message" << endl;
  exit(-1);
}


int main(int argc, char** argv)
{
  int npix = 10000;
  string rawflat_name;
  string scatterflat_name;
  string photflat_name;
  string imgname;
  string ntname;
  
  char c;
  while( (c=getopt(argc, argv, "r:s:p:o:n:")) != -1 ) 
    switch(c)
      {
      case 'r':
	rawflat_name = optarg;
	break;
      case 's':
	scatterflat_name = optarg;
	break;
      case 'p':
	photflat_name = optarg;
	break;
      case 'o':
	imgname = optarg;
	break;
      case 'n':
	ntname = optarg;
	break;
      default:
	usage();
      }
  
  if( rawflat_name == "" ||
      scatterflat_name == "" ||
      photflat_name == "" )
    usage();
  
  FitsImage rawflat(rawflat_name);
  FitsImage scatterflat(scatterflat_name);
  FitsImage photflat(photflat_name);
  
  int nx=rawflat.Nx(), ny=rawflat.Ny();
  assert( (rawflat.Nx() == scatterflat.Nx()) && 
	  (rawflat.Nx() == photflat.Nx()) );
  assert( (rawflat.Ny() == scatterflat.Ny()) && 
	  (rawflat.Ny() == photflat.Ny()) );
  
  int i, sz = nx*ny;
  int dpix = sz / npix;;
  Pixel* rawpix = rawflat.begin();
  Pixel* scatterpix = scatterflat.begin();
  Pixel* photpix = photflat.begin();
  double* vals = new double[npix+1];
  ofstream* ofs = 0;
  if(ntname != "")
    {
      ofs = new ofstream(ntname.c_str());
      (*ofs) << setprecision(12) 
	     << "# i : " << endl
	     << "# val : " << endl
	     << "# end" << endl;
    }
  
  int count = 0;
  for(i=0;i<sz;i+=dpix)
    {
      double r = rawpix[i];
      double s = scatterpix[i];
      double p = photpix[i];
      if( fabs(r)<1.E-6 || fabs(s)<1.E-6 || fabs(p)<1.E-6 ) 
	continue;
      double norm = p / (r * s) ;
      if(count>=npix) continue;
      vals[count++] = norm;
      if(ofs) 
	(*ofs) << i << " "
	       << norm << endl;
    }
  
  double mean, median, sigma;
  Dmean_median_sigma(vals, count, mean, median, sigma);
  cout << " (*) " << mean << " (" << sigma << ") "
       << " [" << median << "]" << endl;
  delete[] vals;
  
  if(imgname != "")
    {
      cout << " (*) imgname=" << imgname << endl;
      FitsImage res(imgname, nx, ny);
      for(i=0;i<sz;i++)
	{
	  double r = rawpix[i];
	  double s = scatterpix[i];
	  double p = photpix[i];
	  if( fabs(r)<1.E-6 || fabs(s)<1.E-6 || fabs(p)<1.E-6 ) 
	    res.begin()[i] = 0.;
	  double norm = p / (r*s);
	  if(norm > 10. || norm<1.E-3) 
	    res.begin()[i] = 0.;
	  else
	    res.begin()[i] = norm;
	}
    }
  
}
