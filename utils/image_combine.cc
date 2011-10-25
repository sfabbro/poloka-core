#include <iostream>

#include "allreducedimage.h"
#include "imagesum.h"


static void usage(const char *progName) {
  cerr << "Usage: " << progName << " [OPTION]... DBIMAGE...\n"
       << "Combine pixels of DBIMAGE\n\n"
       << "    -m METHOD: specify the combining method:\n"
       << "            1  Weighted average (default)\n"
       << "            2  Clipped weighted average\n"
       << "            3  Median\n"
       << "    -w METHOD: specify the weighting method:\n"
       << "             1  Point source optimal (seeing weighted)\n"
       << "             2  Extended source optimal (default)\n"
       << "             3  No global weighting, but using weight maps\n"
       << "             4  No weights at all\n"
       << "    -o DBIMAGE: output name of the combined dbimage (default is 'stack')\n"
       << "    -r DBIMAGE: photometric reference (default is first DBIMAGE)\n";
  exit(EXIT_FAILURE);
}

int main(int nargs, char **args) {

  if (nargs < 2) usage(args[0]);
  if (nargs == 2) { 
    cerr << args[0] << ": need at least 2 images \n";
    return EXIT_FAILURE;
  }

  StackingMethod stackMethod = WeightedAverage;
  WeightingMethod weightMethod = ExtendedSourceOptimal;
  string outName("stack"), phoName("");
  ReducedImageList imList;

  for (int i=1; i<nargs; ++i) {
    char *arg = args[i];
    if (arg[0] != '-') {
      ReducedImageRef im = ReducedImageNew(arg);
      if (!im || !im->IsValid()) { 
	cerr << " not a valid dbimage: " << arg << endl;
	continue;
      }
      imList.push_back(im);
      continue;
    }

    switch (arg[1]) {
    case 'm': { stackMethod = (StackingMethod) atoi(args[++i]); continue; }
    case 'w': { weightMethod = (WeightingMethod) atoi(args[++i]); continue; }
    case 'o': { outName = args[++i]; continue; }
    case 'r': { phoName = args[++i]; continue; }
    default: usage(args[0]);
    }
  }
  if (imList.empty()) { cerr << " no images to combine\n";  return EXIT_FAILURE; }
  if (phoName.empty()) phoName = imList.front()->Name();

  ReducedImage phoRef(phoName);

  cout << " Will combine " << imList.size() << " images:\n"
       << "   - output image    : " << outName << endl
       << "   - photometric ref : " << phoName << endl
       << "   - combine method  : " << name_of_stackingMethod(stackMethod) << endl
       << "   - weight method   : " << name_of_weightingMethod(weightMethod) << endl;

  ImageSum stack(outName, imList, &phoRef, weightMethod, stackMethod);
  stack.Execute(ToTransform(phoRef));

  return EXIT_SUCCESS;
}
