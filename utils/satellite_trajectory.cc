// -*- C++ -*- 
// 
// satellite_trajectory.cc
// 
// 
#include <getopt.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

#include <frame.h>
#include <image.h>
#include <fitsimage.h>



using namespace std;


void usage()
{
  cout << "usage: satellite_trajectory [OPTIONS] <image>" << endl;
  cout << "OPTIONS: " << endl;
  cout << "   -s x,y          specify the starting point" << endl;
  cout << "   -d delta        specify the window range" << endl;
  cout << "   -c <check.fits> check image name" << endl;
  cout << "   -o <traj.list>  satellite track" << endl;
  cout << "   -h              print this message" << endl;
  exit(1);
}



int main(int argc, char** argv)
{
  double x0=0., y0=0.;
  double delta = 5.;
  string im_name, check_im_name, track_nt_name="track.list";
  
  char c;
  while( (c=getopt(argc, argv, "s:d:c:o:h")) != -1 )
    switch(c) 
      {
      case 's':
	{
	  string tmp;
	  istringstream iss(optarg);
	  getline(iss, tmp, ',');
	  x0 = atof(tmp.c_str());
	  getline(iss, tmp);
	  y0 = atof(tmp.c_str());
	}
	break;
      case 'd':
	delta = atof(optarg);
	break;
      case 'c':
	check_im_name = optarg;
	break;
      case 'o':
	track_nt_name = optarg;
	break;
      case 'h':
      default:
	usage();
      }
  if(optind >= argc)
    usage();
  im_name = argv[optind];
  
  double x=x0, y=y0;
  
  FitsImage image(im_name);
  FitsImage* check;
  if(check_im_name != "")
    check = new FitsImage(check_im_name, (Image&)image);
  
  ofstream ofs(track_nt_name.c_str());
  ofs << "# x : " << endl
      << "# y : " << endl
      << "# sky : " << endl
      << "# np : " << endl
      << "# sum : " << endl
      << "# sumx : " << endl
      << "# end" << endl;
  int ix = round(x0);
  
  for(int iy=0; iy<image.Ny(); iy++)
    {
      Pixel mean, sigma;
      double xmin = max(ix-100, 0);
      double xmax = min(ix+100, 1024);
      double ymin = max(iy-100, 0);
      double ymax = min(iy+100, 4096);
      Frame frame(xmin, ymin, xmax, ymax);
      image.SkyLevel(frame, &mean, &sigma);
      
      int np = 0;
      double sum = 0, sumx=0;
      for(ix=x-delta; ix<x+delta; ix++)
	{
	  if(ix<0 || ix>1024) continue;
	  double flx = image(ix,iy)-mean;
	  sumx += flx * ix;
	  sum  += flx;
	  np += 1;
	}
      x = sumx / sum;
      ix = round(x);
      
      if(sum == 0.) break;
      if(check)
	(*check)(ix,iy) = -10;
      
      ofs << x << " " 
	  << iy << " " 
	  << mean << " "
	  << np << " "
	  << sum << " "
	  << sumx << " "
	  << endl;
    }
  
}

