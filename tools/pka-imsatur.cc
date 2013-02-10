#include <cstdio>
#include <poloka/reducedimage.h>

static void usage(const char* progname) {
  cerr << "Usage: " << progname << " [OPTION]...DBIMAGE...\n"
       << "Compute saturation level and create a saturation map\n\n"
       << "    -s : compute a saturation on each amplifier\n"
       << "    -n : print saturation, do not update anything\n"
       << "    -o : overwrite\n\n";
  exit(EXIT_FAILURE);
}

int main(int nargs, char ** args) {

  if (nargs<2) usage(args[0]);

  bool overwrite = false;
  bool sepamps = false;
  bool ok = true;
  bool update = true;

  ReducedImageList imList;

  for (int i=1; i<nargs; ++i)  {
    char *arg = args[i];
    if ((arg[0] != '-')) {
      ReducedImage* im = new ReducedImage(arg);
      if (!im || !im->IsValid()) { 
	cerr << args[0] << ": " << arg << "is not a valid dbimage\n";
	ok = false;
	continue;
      }
      imList.push_back(im);
      continue;
    }
    switch (arg[1]) {
    case 'o': overwrite = true; break;
    case 'n': update = false; break;
    case 's': sepamps = true; break;
    default: usage(args[0]);
    }
  }
  
  for(ReducedImageIterator it = imList.begin(); it != imList.end(); ++it) {
    ReducedImageRef current = *it;
    if (update) {
      if (overwrite && current->HasSatur())
	ok = remove((current->FitsSaturName()).c_str());
      ok = current->MakeSatur();
    } else {
      
    }    
  }

  return ok? EXIT_SUCCESS : EXIT_FAILURE;
}
