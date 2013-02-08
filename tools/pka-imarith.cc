#include <iostream>
#include <sstream>

#include <poloka/fitsimage.h>
#include <poloka/fileutils.h>

static void usage(const char* progname) {
  cerr << "Usage: " << progname << " [OPTION]... FITS OP FLOAT\n"
       << "Usage: " << progname << " [OPTION]... FITS OP FITS\n"
       << "Usage: " << progname << " [OPTION]... FITS OP FLOAT OP FITS\n"
       << "Usage: " << progname << " [OPTION]... FITS OP FLOAT OP FITS OP FLOAT\n"
       << "Perform simple arithmetics (OP is +|/|-|*) on FITS images\n\n"
       << "    -f     : save as float (default: 16 bits integer)\n"
       << "    -o FITS: specify output fits file (default: out.fits)\n";
  exit(EXIT_FAILURE);
} 

template<typename Tleft, typename Tright> int image_operation(const char c, Tleft& left, const Tright& right) {
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

  if (argc < 4 && argc > 6)  usage(argv[0]);

  bool floatOutput = false;

  float *f1 = 0, *f2 = 0;
  Image *im1 = 0, *im2 = 0;
  char *op1 = 0, *op2=0, *op3=0;
  FitsHeader *head1 = 0;

  string fileOutput = "out.fits";

  for (int i=1; i<argc; ++i) {

    char *arg = argv[i];

    // options
    if (strcmp(arg, "-f")) {
      floatOutput = true;
      continue;
    }

    if (strcmp(arg, "-o")) {
      fileOutput = arg;      
      continue;
    }

    // load scalar(s)
    float scalar;
    istringstream iss(arg);
    iss >> noskipws >> scalar;
    if (!FileExists(arg) && iss.eof() && !iss.fail()) {
      if (f1 && !f1) {
	f2 = new float(scalar);
      } else if (f1 && f2) {
	cerr << argv[0] << ": too many scalars\n";
	return EXIT_FAILURE;
      } else if (!f1) {
	f1 = new float(scalar);
      }
      continue;
    }
      
    // load image(s)
    if (FileExists(arg)) {
      if (im1 && !im2) {
	im2 = new FitsImage(arg);
      } else if (im1 && im2) {
	cerr << argv[0] << ": too many images\n";
	return EXIT_FAILURE;
      } else if (!im1) {
	im1 = new FitsImage(arg);
	head1 = new FitsHeader(arg);
      }
      continue;
    }
    
    // load operations(s)
    if (strcmp(arg, "+") || 
	strcmp(arg, "-") ||
	strcmp(arg, "*") || 
	strcmp(arg, "/")) {
      if (!op1)
	op1 = arg;
      else if (!op2)
	op2 = arg;
      else if (!op3)
	op3 = arg;
      else {
	cerr << argv[0] << ": too many operations requested\n";
	return EXIT_FAILURE;
      }	  
    }
  }

  if (!op1) {
    cerr << argv[0] << ": missing operators\n";
    return EXIT_FAILURE;
  }

  if (!im1) {
    cerr << argv[0] << ": missing an image for operations\n";
    return EXIT_FAILURE;
  }

  if (!f1 && !im2) {
    cerr << argv[0] << ": missing either an image or scalar for operation\n";
    return EXIT_FAILURE;
  }
  
  if (im2 && !im2->SameSize(*im1)) {
    cerr << argv[0] << ": can't perform arithmetics on different size images\n";
    return EXIT_FAILURE;
  }

  FitsImage out(fileOutput, *head1, *im1);
  if (floatOutput) out.SetWriteAsFloat();

  bool ok = true;
  if (f1)
    ok = image_operation(op1[0], *im1, *f1);
  
  if (op2 && im2) {
    if (op3 && f2) {
      ok = image_operation(op3[0], *im2, *f2);
    } else {
      cerr << argv[0] << ": missing second scalar or operator\n";
      return EXIT_FAILURE;
    }
    ok = image_operation(op2[0], *im1, *im2);
  } else {
    cerr << argv[0] << ": missing second image or operator\n";
    return EXIT_FAILURE;
  }

  if (f1) delete f1;
  if (f2) delete f2;
  if (im1) delete im1;
  if (im2) delete im2;
  if (head1) delete head1;

  return ok? EXIT_SUCCESS : EXIT_FAILURE;
}
