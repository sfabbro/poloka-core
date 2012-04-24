#include <cmath>
#include <cstdlib>
#include <ctime>

#include "allreducedimage.h"
#include "basestar.h"
#include "imagepsf.h"
#include "fileutils.h"
#include "fitsimage.h"

static void usage(const char* prog) {
  cerr << prog << " <file with ra dec mag> <dbimage>...<dbimage>\n";
  exit(EXIT_FAILURE);
}

void random_init(int& seed) {
  if (seed == 0)
    seed = time(NULL);
  srand(seed);
}


static double random_uniform() {
  return double(rand()) / RAND_MAX;
}


static double gammln(const double& xx) {
  static const double cof[6]= 
    { 76.18009173, -86.50532033, 24.01409822,
      -1.231739516, 0.120858003e-2, -0.536382e-5};

  double x = xx - 1.0;
  double tmp = x + 5.5;
  tmp -= (x+0.5)*log(tmp);
  double ser = 1.0;
  for (int j=0; j<6; j++)
    ser += cof[j]/(x+=1.0);
  return log(2.50662827465*ser) - tmp;
}

static double random_poisson(const double& xm) {
  static double sq, alxm, g, oldm=-1.0;
  double em,t,y;
  if (xm < 12.0)  {
    if (xm != oldm) {
      oldm = xm;
      g = exp(-xm);
    }
    em = -1.0;
    t = 1.0;
    do {
      em += 1.0;
      t *= random_uniform();
    } while (t > g);
  } else {
    if (xm != oldm) {
      oldm = xm;
      sq = sqrt(2.0*xm);
      alxm = log(xm);
      g = xm*alxm - gammln(xm+1.0);
    }
    do {
      do {
        y = tan(M_PI*random_uniform());
        em = sq*y+xm;
      } while (em < 0.0);
      em = floor(em);
      t = 0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
    } while (random_uniform() > t);
  }

  return em;
}

struct ImageAddStar {
  BaseStarList stars;
  ImageAddStar(const char* filename) {
    stars.read(filename);
  }

  void operator () (const string& name) {
    ReducedImageRef rim = ReducedImageNew(name);
    if (!rim->IsValid()) { 
      cerr << " not a valid dbimage: " << name << endl;
      return;
    }

    ImagePSF psf(*rim, false);
    GtransfoRef wcs = rim->RaDecToPixels();
    string outfile = rim->Dir() + "planted.fits";
    CopyFile(rim->FitsName(), outfile);
    FitsImage image(outfile, RW);
    double zp = rim->AnyZeroPoint();
    double gain = rim->Gain();
    double sig2 = rim->SigmaBack();
    sig2 *= sig2;
    cout << " Planting sources with zeropoint = " << zp << endl;
    cout << " # x y flux s/n mag\n";
    for (BaseStarCIterator it=stars.begin(); it != stars.end(); ++it) {
      const BaseStar* s = *it;
      double flux = pow(10,0.4*(zp - s->flux));
      double x,y;
      wcs->apply(s->x, s->y, x, y);
      int istart, jstart, iend, jend;
      psf.StampLimits(x, y, istart, iend, jstart, jend);
      if ((istart<iend) && (jstart<jend)) {
	Image star(iend-istart, jend-jstart);
	double signal = 0;
	double norm = 0;
	double vmin = 0;
	for (int i=istart; i<iend; ++i)
	  for (int j=jstart; j<jend; ++j) {
	    double p = psf.PSFValue(x, y, i, j);
	    double adu = flux * p;
	    signal += adu;
	    double signal_poisson = random_poisson(adu*gain)/gain;
	    norm += signal_poisson;
	    star(i-istart,j-jstart) = signal_poisson;
	    vmin += p*p / (adu/gain + sig2);
	  }
	// make sure total(added signal) = flux after adding poisson noise
	norm = signal / norm;
	for (int i=istart; i<iend; ++i)
	  for (int j=jstart; j<jend; ++j)
	    image(i,j) += norm * star(i-istart,j-jstart);

	cout << x << " " << y
	     << " " << flux
	     << " " << signal*sqrt(vmin)
	     << " " << s->flux << endl;
      }
    }
  }
};

int main( int nargs, char **args){
  if (nargs <=1) usage(args[0]);

  ImageAddStar imAddStar(args[1]);
  list<string> imList;

  for (int i=2; i < nargs; ++i)
    imList.push_back(args[i]);

  int seed=0;
  random_init(seed);
  for_each(imList.begin(), imList.end(), imAddStar);

  return EXIT_SUCCESS;
}
