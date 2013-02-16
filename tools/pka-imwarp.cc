#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <poloka/swarpstack.h>

static void usage(const char* progname) {
  cerr << "Usage: " << progname << " [OPTION]... DBIMAGE...\n" 
       << "Wrapper around SWarp image warping and stacking software\n\n"
       << "     -s         : stack subtractions instead of calibrated images\n"
       << "     -o DBIMAGE : output stack DBIMAGE (default: stack)\n"
       << "     -r DBIMAGE : specify a astrometric/photometric reference\n"
       << "     -c FILE    : specify a SWarp configuration file instead of default\n"
       << "     -i FILE    : file with a list of DBIMAGEs instead of argument list\n\n";
  exit(EXIT_FAILURE);
}

int main(int nargs, char **args) {

  if (nargs <= 1) usage(args[0]);

  ReducedImageList imList;
  string outName = "stack";
  string cardsName;
  ReducedImageRef ref;
  map<string,Point> dcrval;
  DbImageKind imageType = Calibrated;

  for (int i=1; i<nargs; ++i) {
    char *arg = args[i];
    if (arg[0] == '-') {
      switch (arg[1]) {
      case 'r': 
	ref = new ReducedImage(args[++i]); 
	if (ref && !ref->HasImage()) {
	  cerr << args[0] << ": the reference" << ref->Name() << " does not have a fits image\n";
	  return EXIT_FAILURE;
	}
	break;
      case 's': imageType = Subtracted; break;
      case 'o': outName = args[++i]; break;
      case 'c': cardsName = args[++i]; break;
      case 'i': {
	ifstream ifs(args[++i]);
	if (ifs.is_open()) {
	  while (ifs.good()) {
	    string line;
	    getline(ifs, line);
	    if (line.empty() || line[0] == '#') continue;
	    string name;
	    istringstream iline(line);
	    iline >> name;
	    ReducedImageRef ri = new ReducedImage(name);
	    if (ri && ri->IsValid()) 
	      imList.push_back(ri);
	    else 
	      cerr << args[0] << ": " << arg << " is not a valid dbimage\n";
	    Point pt;
	    if (iline >> pt.x >> pt.y)
	      dcrval[name] = pt;
	  }
	}
      }
	break;
      default: usage(args[0]);
      }
      continue;
    } else {
      ReducedImageRef ri = new ReducedImage(args[i]);
      if (ri && ri->IsValid()) 
	imList.push_back(ri);
      else 
	cerr << args[0] << ": " << arg << " is not a valid dbimage\n";
    }
  }
  
  if (imList.empty()) {
    cerr << args[0] << ": no input images provided\n";
    return EXIT_FAILURE;
  }
  
  if (outName.empty()) {
    cerr << args[0] << ": no output name provided\n";
    return EXIT_FAILURE;
  }
  
  if (!cardsName.empty())
    SetSwarpCardsName(cardsName);
  
  Frame frame;
  if (ref) frame = ref->PhysicalSize();

  SwarpStack ss(outName, imList, ref, frame);
  bool ok = true;
  ss.SetSwarpType(imageType);
  ss.dcrval = dcrval;
  ok = ss.MakeFits();
  ok = ss.MakeSatur();
  ok = ss.MakeCatalog();

  return ok? EXIT_SUCCESS : EXIT_FAILURE;
}
