#include <cstdio>

#include <poloka/reducedimage.h>
#include <poloka/polokaexception.h>

static void usage(char * progname) {
  cerr << "Usage: " << progname << " DBIMAGE...\n"
       << "Produce a binary map of the satellites\n\n"
       << "   -o: overwrite\n";
  exit(EXIT_FAILURE);
}

int main(int nargs, char ** argv) {
  if(nargs ==1) usage(argv[0]);

  bool overwrite = false;
  ReducedImageList imList;
  bool ok = true;

  for (int i=1; i< nargs; i++) {
    char *arg = argv[i];
    if (arg[0] != '-') {
      ReducedImage *current = new ReducedImage(argv[i]);
      if (current && current->IsValid())
	imList.push_back(current);
      else {
	cerr << argv[0] << ": " << argv[i] << " not found\n";
	ok = false;
      }
      continue;
    }
    switch (arg[1]) {
    case 'o': overwrite = true; break;
    default:  usage(argv[0]);
    }
  }

  for (ReducedImageIterator it = imList.begin(); it != imList.end(); ++it) {
    try {
      ReducedImage* current = *it;
      if (overwrite && current->HasSatellite())
	ok = remove((current->FitsSatelliteName()).c_str());
      ok = current->MakeSatellite();
    } catch(PolokaException p) {
      p.PrintMessage(cerr);
      ok = false;
    }
  }
  
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
