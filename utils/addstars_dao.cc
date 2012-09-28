#include <cmath>
#include <cstdlib>
#include <ctime>

#include "allreducedimage.h"
#include "reducedimage.h"
#include "daostar.h"
#include "daophotio.h"

#include "fileutils.h"
#include "fitsimage.h"

static void usage(const char* prog) {
  cerr << prog << " <file with id ra dec mag> <dbimage>...<dbimage>\n";
  exit(EXIT_FAILURE);
}

struct genstar {
  string id;
  double ra,dec,mag;
};

static void read_generated(const char* filename, list<genstar>& stars) {
  ifstream in(filename);
  char c;
  while (in >> c) {
    in.unget();
    genstar star;
    in >> star.id >> star.ra >> star.dec >> star.mag;
    stars.push_back(star);    
  }
}


struct ImageAddStar {
  list<genstar> stars;

  ImageAddStar(const char* filename) {
    read_generated(filename, stars);
  }

  void operator () (const string& name) {
    ReducedImageRef rim = ReducedImageNew(name);
    if (!rim->IsValid()) { 
      cerr << " not a valid dbimage: " << name << endl;
      return;
    }

    GtransfoRef wcs = rim->RaDecToPixels();
    if (!wcs) {
      cerr << " error in converting ra dec to x y\n";
      return;
    }
    double zp = rim->AnyZeroPoint();
    Frame frame = rim->UsablePart();
    cout << " # Planting sources with zeropoint = " << zp << endl;
    cout << " # x y flux mag\n";
    DaoStarList daostars;
    int i = 1;
    for (list<genstar>::iterator it=stars.begin(); it != stars.end(); ++it) {
      double flux = pow(10, 0.4*(zp - it->mag));
      double x,y;
      wcs->apply(it->ra, it->dec, x, y);
      if (!frame.InFrame(x,y)) continue;
      DaoStar *daostar = new DaoStar;
      daostar->num = i++;
      daostar->x = x;
      daostar->y = y;
      daostar->flux = flux;
      daostar->sky = 0;      
      daostars.push_back(daostar);
      cout << x << " " << y << " " << flux << " " << it->mag << endl;
    }
    WriteDaoList(*rim, "planted.lst", daostars);
  }
};

int main( int nargs, char **args){
  if (nargs <=1) usage(args[0]);

  ImageAddStar imAddStar(args[1]);
  list<string> imList;

  for (int i=2; i < nargs; ++i)
    imList.push_back(args[i]);

  for_each(imList.begin(), imList.end(), imAddStar);

  return EXIT_SUCCESS;
}
