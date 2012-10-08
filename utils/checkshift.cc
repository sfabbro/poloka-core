#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "reducedimage.h"

static void usage(const char* prog) {
  cerr << prog << " REFDB SHIFTFILE\n";
  exit(EXIT_FAILURE);
}

struct genstar {
  string id;
  double ra,dec,mag;
  bool operator == (const string &right) const { return id == right; }
};

static double sq(const double& x) { return x*x; }

struct MatchStar {

  list<genstar> refStars;
  double mjd;

  MatchStar(const string& name) {

    ReducedImage ref(name);
    if (!ref.IsValid()) { 
      cerr << " not a valid dbimage: " << name << endl;
      return;
    }    

    ifstream in((ref.Dir()+"/planted.list").c_str());
    if (!in) {
      cerr << " error reading " << ref.Dir()+"/planted.list \n";
      return;
    }
    char c;
    while (in >> c) {
      in.unget();
      genstar star;
      in >> star.id >> star.ra >> star.dec >> star.mag;
      refStars.push_back(star);
    }
  }

  void operator () (const string& name) {

    ReducedImage im(name);
    if (!im.IsValid()) { 
      cerr << " not a valid dbimage: " << name << endl;
      return;
    }

    GtransfoRef wcs = im.RaDecToPixels();
    ifstream in((im.Dir()+"/planted.list").c_str());
    GtransfoRef wcsinv = im.PixelsToRaDec();
    //ofstream out((im.Name()+".shift").c_str());
    //double mdist=0, sdist=0, mdx=0, sdx=0, mdy=0, sdy=0;
    //double ra, dec, mag, x, y, xref, yref, rainv, decinv;
    double ra, dec, mag, x, y, xref, yref;
    double mdx=0, mdy=0;
    Frame frame = im.UsablePart();
    char c;
    double nstars = 0;
    while (in >> c) {
      in.unget();
      string id;
      in >> id >> ra >> dec >> mag;
      wcs->apply(ra, dec, x, y);
      if (!frame.InFrame(x,y)) continue;
      list<genstar>::iterator rit = find(refStars.begin(), refStars.end(), id);
      if (rit != refStars.end()) {
	wcs->apply(rit->ra, rit->dec, xref, yref);
	//wcsinv->apply(x+dx, y+dy, rainv, decinv);
	//out << setiosflags(ios::fixed)
	//    << setw(4)
	//    << id << " "
	//    << setprecision(6)
	//    << rit->ra << " "
	//    << rit->dec << " "
	//    << rainv << " " 
	//    << decinv << " "
	//    << setprecision(2)
	//    << mag << " "
	//    << setprecision(3)
	//    << setw(9)
	//    << xref << " "
	//    << setw(9)
	//    << yref << " "
	//    << setw(9)
	//    << x << " "
	//    << setw(9)
	//    << y
	//    << endl;
	//double dist = sqrt( sq((ra-rainv)*cos(M_PI*decinv/180.)) + sq(dec-decinv)) * 3600.;
	//mdist += dist;
	//sdist += sq(dist);
	mdx += xref - x;
	mdy += yref - y;
	//sdx += sq(x-xref);
	//sdy += sq(y-yref);
	nstars++;
      }
    }

    //mdist /= nstars;
    //sdist = sdist/nstars - sq(mdist);
    mdx /= nstars;
    //sdx = sdx/nstars - sq(mdx);
    mdy /= nstars;
    //sdy = sdy/nstars - sq(mdy);
    
    cout << name << " " 
	 << setiosflags(ios::fixed)
	 << setprecision(5) 
	 << setw(15)
	 << mjd << " "
	 << setprecision(3)
	 << setw(9)
      // << mdist 
      //<< "(" << sdist << ") "
	 << setw(9)
	 << mdx << " " 
      // << "(" << sdx << ") "
	 << setw(9)
	 << mdy << " "
      // << "(" << sdy << ") "
	 << endl;

    in.close();
    in.open((im.Dir()+"/planted.list").c_str());
    double mdist = 0, sdist = 0;
    while (in >> c) {
      in.unget();
      string id;
      in >> id >> ra >> dec >> mag;
      wcs->apply(ra, dec, x, y);
      if (!frame.InFrame(x,y)) continue;
      list<genstar>::iterator rit = find(refStars.begin(), refStars.end(), id);
      if (rit != refStars.end()) {
	wcsinv->apply(x+mdx, y+mdy, ra, dec);
	double dist = sqrt( sq((ra-rit->ra)*cos(M_PI*dec/180.)) + sq(dec-rit->dec)) * 3600./0.185;
	mdist += dist;
	sdist += sq(dist);
      }
    }
    cerr << " mean distance error " << fixed 
	 << mdist/nstars << " (" << sdist/nstars - sq(mdist/nstars) << ")\n";
  }  
};

int main(int nargs, char **args){
  if (nargs != 3) usage(args[0]);

  MatchStar matchStar(args[1]);
  if (matchStar.refStars.empty()) return EXIT_FAILURE;
  ifstream in(args[2]);
  char c;
  while (in >> c) {
    in.unget();
    string name;
    double mjd, dx, dy;
    in >> name >> matchStar.mjd >> dx >> dy;
    matchStar(name);
  }
  return EXIT_SUCCESS;
}
