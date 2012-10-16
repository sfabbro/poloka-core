#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "reducedimage.h"

static void usage(const char* prog) {
  cerr << prog << " REFDB DBIMAGE...\n";
  exit(EXIT_FAILURE);
}

struct genstar {
  string id;
  double ra,dec,mag;
  bool operator == (const string &right) const { return id == right; }
};

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
	mdx += xref - x;
	mdy += yref - y;
	nstars++;
      }
    }
    mdx /= nstars;
    mdy /= nstars;
    
    cout << name << " " 
	 << setiosflags(ios::fixed)
	 << setprecision(5) 
	 << setw(15)
	 << mjd << " "
	 << setprecision(3)
	 << setw(9)
	 << setw(9)
	 << mdx << " " 
	 << setw(9)
	 << mdy << " "
	 << endl;
  }  
};

int main(int nargs, char **args) {

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
