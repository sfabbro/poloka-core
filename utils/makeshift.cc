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

static double compute_step(const ReducedImageList& imlist) {
  size_t n = imlist.size();
  vector<double> seeing(n);
  size_t i = 0;
  for (ReducedImageCIterator it = imlist.begin(); it != imlist.end(); ++it, ++i)
    seeing[i] = (*it)->Seeing() * (*it)->PixelSize() / 3600;
  // median seeing is resulting step size
  sort(seeing.begin(), seeing.end());
  return (n & 1) ? seeing[n/2] : (seeing[n/2-1] + seeing[n/2])*0.5;  
}

struct MatchStar {

  list<genstar> refStars;
  string refname;
  double step;
  MatchStar(const string& name) : refname(name) {

    ReducedImage ref(name);
    if (!ref.IsValid()) { 
      cerr << " MatchStar: " << name << " is not a valid dbimage\n";
      return;
    }    

    ifstream in((ref.Dir()+"/planted.list").c_str());
    if (!in) {
      cerr << " MatchStar: error reading " << ref.Dir() + "/planted.list \n";
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
      cerr << " MatchStar: " << name << " is not a valid dbimage\n";
      return;
    }

    GtransfoRef wcs = im.RaDecToPixels();
    ifstream in((im.Dir()+"/planted.list").c_str());
    if (!in) {
      cerr << " MatchStar: error reading " << im.Dir() + "/planted.list \n";
      return;
    }
    Frame frame = im.UsablePart();
    double minrashift = 1e10, mindecshift = 1e10;
    double maxrashift = 1e-10, maxdecshift = -1e10;
    double rashift, decshift;

    char c;
    double nstars = 0;
    while (in >> c) {
      in.unget();
      string id;
      double ra, dec, mag, x, y;
      in >> id >> ra >> dec >> mag;
      wcs->apply(ra, dec, x, y);
      if (!frame.InFrame(x,y)) continue;
      list<genstar>::iterator rit = find(refStars.begin(), refStars.end(), id);
      if (rit != refStars.end()) {
	rashift = rit->ra - ra;
	decshift= rit->dec - dec;
	if (rashift < minrashift) minrashift = rashift;
	else if (rashift > maxrashift) maxrashift = rashift;
	if (decshift < mindecshift) mindecshift = decshift;
	else if (decshift > maxdecshift) maxdecshift = decshift;
	nstars++;
      }
    }
    in.close();
    if (nstars == 0) return;

    rashift = minrashift;
    decshift = mindecshift;
    while (rashift <= maxrashift) {
      rashift += step;
      while (decshift <= maxdecshift) {
	decshift += step;
	string filename = "shifts_" + grid.str() + ".dat";
	ofstream out(filename().c_str(), ios::app);
	out << name << " "
	    << setiosflags(ios::fixed)
	    << setprecision(10)
	    << setw(9)
	    << rashift << " "
	    << decshift << " "
	    << endl;
      }
    }
  }
};

int main(int nargs, char **args) {

  if (nargs < 2) usage(args[0]);
  ReducedImageList imList;

  MatchStar matchStar(args[1]);
  
  if (matchStar.refStars.empty())
    return EXIT_FAILURE;

  for (int i=2; i<nargs; ++i)
    imList.push_back(args[i]);

  matchStar.step = compute_step(imList);
  foreach(imList.begin(), imList.end(), matchStar);
  
  return EXIT_SUCCESS;
}
