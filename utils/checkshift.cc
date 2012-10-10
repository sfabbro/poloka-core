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
  double ra,dec,mag,x,y;
  bool operator == (const string &right) const { return id == right; }
};

static double sq(const double& x) { return x*x; }

struct MatchStar {

  list<genstar> refStars;
  double mjd;
  GtransfoRef ref2xy;

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
    ref2xy = ref.RaDecToPixels();
    while (in >> c) {
      in.unget();
      genstar star;
      in >> star.id >> star.ra >> star.dec >> star.mag;
      ref2xy->apply(star.ra, star.dec, star.x, star.y);
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
    double mdx=0, sdx=0, mdy=0, sdy=0;
    double ra, dec, mag, x, y, xref, yref;
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
	mdx += xref - x;
	mdy += yref - y;
	sdx += sq(x-xref);
	sdy += sq(y-yref);
	nstars++;
      }
    }

    mdx /= nstars;
    sdx = sqrt(sdx/nstars - sq(mdx));
    mdy /= nstars;
    sdy = sqrt(sdy/nstars - sq(mdy));
    
    cout << name << " " 
	 << setiosflags(ios::fixed)
	 << setprecision(5) 
	 << setw(15)
	 << mjd << " "
	 << setprecision(3)
	 << setw(9)
	 << mdx << " " 
	 << "(" << sdx << ") "
	 << setw(9)
	 << mdy << " "
	 << "(" << sdy << ") ";

    in.close();

    in.open((im.Dir()+"/planted.list").c_str());
    GtransfoRef wcsinv = im.PixelsToRaDec();
    ofstream out((im.Name()+".shift").c_str());
    mdx = mdy = sdx = sdy = 0;
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
	ref2xy->apply(ra, dec, x, y);
	double dist = sqrt( sq((ra-rit->ra)*cos(M_PI*dec/180.)) + sq(dec-rit->dec)) * 3600./0.186;
	double dx = rit->x - x;
	double dy = rit->y - y;
	out << setiosflags(ios::fixed)
	    << setprecision(2) 
	    << rit->mag << " " 
	    << setprecision(6) 
	    << setw(12)
	    << rit->ra << " "
	    << rit->dec << " "
	    << setprecision(6) 
	    << setw(10)
	    << dx << " "
	    << setprecision(6) 
	    << setw(10)
	    << dy << " "
	    << setprecision(6) 
	    << setw(10)
	    << dist << endl;
	mdx += dx;
	sdx += sq(dx);
	mdy += dy;
	sdy += sq(dy);
	mdist += dist;
	sdist += sq(dist);
      }
    }

    mdx /= nstars;
    sdx = sqrt(sdx/nstars - sq(mdx));
    mdy /= nstars;
    sdy = sqrt(sdy/nstars - sq(mdy));
    mdist /= nstars;
    sdist = sqrt(sdist/nstars - sq(mdist));

    cout << mdx << " (" << sdx << ") "
	 << mdy << " (" << sdy << ") "
	 << mdist << " (" << sdist << ")\n";
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
    string name; double dum;
    in >> name >> matchStar.mjd >> dum >> dum;
    matchStar(name);
  }
  return EXIT_SUCCESS;
}
