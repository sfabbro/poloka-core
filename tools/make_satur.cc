#include "allreducedimage.h"
#include "reducedimage.h"
#include <cstdio>

static void usage(const char* progname) {
  cerr << "Usage: " << progname << " [OPTION]...DBIMAGE...\n"
       << "Compute saturation level and create a saturation map\n\n"
       << "    -o : overwrite\n";
  exit(EXIT_FAILURE);
}

int main(int nargs, char ** args) {

  if (nargs<2) usage(args[0]);

  bool overwrite = false;
  ReducedImageList imList;

  for (int i=1; i<nargs; ++i)  {
    char *arg = args[i];
    if ((arg[0] != '-')) {
      ReducedImage* im = new ReducedImage(arg);
      if (!im || !im->IsValid()) { 
	cerr << arg << ": not a valid dbimage\n";
	continue;
      }
      imList.push_back(im);
      continue;
    }
    switch (arg[1]) {
    case 'o': overwrite = true; break;
    default: usage(args[0]);
    }
  }
  
  for(ReducedImageIterator it = imList.begin(); it != imList.end(); ++it) {
    ReducedImageRef current = *it;
    if (overwrite && current->HasSatur())
      remove((current->FitsSaturName()).c_str());
    current->MakeSatur();
  }

  return EXIT_SUCCESS;
}
