#include <iostream>
#include <sstream>

#include <poloka/fitsimage.h>
#include <poloka/fileutils.h>

static void usage(const char* progname) {
  cerr << "Usage: " << progname << " [OPTION]... FITS|FLOAT OP FLOAT|FITS\n"
       << "Perform arithmetics (OP is +|/|-|*) on FITS images with FLOAT scalar\n\n"
       << "    -f     : save as float (default: 16 bits integer)\n"
       << "    -o FITS: specify output fits file (default: out.fits)\n\n";
  exit(EXIT_FAILURE);
} 

template<typename Tright> int image_operation(const char c, Image& left, const Tright& right) {
  switch (c) {
  case '+':
    left += right;
    break;
  case '*':
    left *= right;
    break;
  case '/':
    left /= right;
    break;
  case '-':
    left -= right;
    break;
  default:
    return false;
  }
  return true;
}

int main(int argc, char** argv) {

  if (argc < 4 || argc > 6)  usage(argv[0]);

  bool floatOutput = false;

  float *f = 0;
  Image *im1 = 0;
  Image *im2 = 0;
  char op = '0';
  FitsHeader *head1 = 0;
  bool imagefirst = false;

  string fileOutput = "out.fits";

  for (int i=1; i<argc; ++i) {

    char *arg = argv[i];

    // options
    if (strcmp(arg, "-f") == 0) {
      floatOutput = true;
      continue;
    }

    if (strcmp(arg, "-o") == 0) {
      fileOutput = arg;      
      continue;
    }

    // load scalar
    float scalar;
    istringstream iss(arg);
    iss >> noskipws >> scalar;
    if (!FileExists(arg) && iss.eof() && !iss.fail()) {
      if (!f) {
	f = new float(scalar);
	if (im1) imagefirst = true;
      } else {
	cerr << argv[0] << ": too many scalars\n";
	return EXIT_FAILURE;
      }
      continue;
    }
      
    // load image(s)
    if (FileExists(arg)) {
      if (im1 && !im2) {
	im2 = new FitsImage(arg);
      } else if (!im1) {
	im1 = new FitsImage(arg);
	head1 = new FitsHeader(arg);
      } else if (im1 && im2) {
	cerr << argv[0] << ": too many images\n";
	return EXIT_FAILURE;
      }
      continue;
    }
    
    // load operations(s)
    if (strcmp(arg, "+") == 0 || 
	strcmp(arg, "-") == 0 ||
	strcmp(arg, "*") == 0 || 
	strcmp(arg, "/") == 0) {
      if (op == '0')
	op = arg[0];
      else {
	cerr << argv[0] << ": too many operations\n";
	return EXIT_FAILURE;
      }
    }
  }

  if (op == '0') {
    cerr << argv[0] << ": missing operator\n";
    return EXIT_FAILURE;
  }

  if (!im1) {
    cerr << argv[0] << ": missing an image for operations\n";
    return EXIT_FAILURE;
  }

  if (!f && !im2) {
    cerr << argv[0] << ": missing either an image or scalar for operation\n";
    return EXIT_FAILURE;
  }
  
  if (im2 && !im2->SameSize(*im1)) {
    cerr << argv[0] << ": can't perform arithmetics on images with different size\n";
    return EXIT_FAILURE;
  }

  bool ok = true;
  if (im1 && im2)
    ok = image_operation(op, *im1, *im2);
  else if (im1 && f) {
    if (!imagefirst) {
      if (op == '-') {
	*im1 *= -1.;
	op = '+';
      } else if (op == '/') {
	*im1 = 1./(*im1);
	op = '*';
      }
    }
    ok = image_operation(op, *im1, *f);
  } else {
    cerr << argv[0] << ": wrong arguments\n";
    return EXIT_FAILURE;
  }
  
  if (f) delete f;
  if (im2) delete im2;

  FitsImage out(fileOutput, *head1, *im1);
  if (floatOutput) out.SetWriteAsFloat();

  if (im1) delete im1;
  if (head1) delete head1;

  return ok? EXIT_SUCCESS : EXIT_FAILURE;
}
