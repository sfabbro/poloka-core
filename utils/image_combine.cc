#include <iostream>

#include "allreducedimage.h"
#include "imagesum.h"


static void usage(const char *progName) {
  cerr << "Usage: " << progName << " [OPTIONS] <DbImage(s)>\n"
       << "  combine pixels (do not resample)\n"
       << "  [OPTIONS] are \n"
       << "     -m <n>: specify the combining method:\n"
       << "              1  Weighted average (default)\n"
       << "              2  Clipped weighted average\n"
       << "              3  Median\n"
       << "     -w <n>: specify the weight method:\n"
       << "              1  Point source optimal\n"
       << "              2  Extended source optimal (default)\n"
       << "              3  No global weighting\n"
       << "              4  No weights at all\n"
       << "     -o <name>: name of the combined dbimage (default is 'stack')\n"
       << "     -r <name>: name of the photometric reference (default is first image)\n";
}

int main(int nargs, char **args) {

  if (nargs < 2) { usage(args[0]); return EXIT_SUCCESS; }
  if (nargs == 2) { 
    cerr << " combine at least 2 images \n";
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
    default: usage(args[0]); return EXIT_FAILURE;
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
