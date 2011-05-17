#include "allreducedimage.h"
#include "reducedimage.h"
#include <cstdio>

static void usage(const string &exec) {
  cerr << "Usage: " << exec << " <DbImages...>\n"
       << "    -o : overwrite\n";
}

int main(int nargs, char ** args) {

  if (nargs<2) { usage(args[0]); return EXIT_FAILURE; }

  bool overwrite = false;
  ReducedImageList imList;

  for (int i=1; i<nargs; ++i)  {
    char *arg = args[i];
    if ((arg[0] != '-')) {
      ReducedImageRef im = ReducedImageNew(arg);
      if (!im || !im->IsValid()) { 
	cerr << " not a valid dbimage: " << arg << endl;
	continue;
      }
      imList.push_back(im);
      continue;
    }
    switch (arg[1]) {
    case 'o': overwrite = true; break;
    default: usage(args[0]); return EXIT_FAILURE;
    }
  }
  
  for(ReducedImageIterator it = imList.begin(); it != imList.end(); ++it) {
    ReducedImage *current = *it;
    if (overwrite && current->HasSatur())
      remove((current->FitsSaturName()).c_str());
    current->MakeSatur();
  }

  return EXIT_SUCCESS;
}
