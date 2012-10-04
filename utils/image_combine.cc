#include <iostream>

#include "allreducedimage.h"
#include "imagesum.h"
#include "polokaexception.h"

static void usage(const char *progName) {
  cerr << "Usage: " << progName << " [OPTION]... DBIMAGE...\n"
       << "Combine pixels of DBIMAGE\n\n"
       << "   -o DBIMAGE: output name of the combined dbimage (default is 'stack')\n"
    //       << "   -b METHOD : background subtraction method:\n"
    //       << "               1  On the fly estimation\n"
    //       << "               2  No background subtraction (default)\n"
       << "   -c METHOD : pixel combining method:\n"
       << "               1  Weighted average (default)\n"
       << "               2  Clipped weighted average\n"
       << "               3  Median\n"
       << "               4  Adaptive weighted average\n"
       << "   -w METHOD : weighting method:\n"
       << "               1  Weight maps & seeing weighted\n"
       << "               2  Weight maps & average S/N weighted (default)\n"
       << "               3  Weight maps\n"
       << "               4  No weighting\n"
       << "   -p METHOD : photometric scaling method\n"
       << "               1  Use zero point ZP_PHOT (default)\n"
       << "               2  Compute star flux ratio with reference\n"
       << "               5  No photometric scaling\n"
       << "   -z THING  : photometric reference: zp, dbimage or catalog (default is 30)\n\n";
  exit(EXIT_FAILURE);
}

int main(int nargs, char **args) {

  if (nargs < 2) usage(args[0]);
  if (nargs == 2) { 
    cerr << args[0] << ": need at least 2 images\n";
    return EXIT_FAILURE;
  }

  StackingMethod stackMethod = WeightedAverage;
  WeightingMethod weightMethod = ExtendedSourceOptimal;
  //BackgroundMethod backMethod = OnTheFly;
  PhotoScalingMethod scaleMethod = TotalLeastSquares;

  string outName("stack"), phoRef("30");
  ReducedImageList imList;

  for (int i=1; i<nargs; ++i) {
    char *arg = args[i];
    if (arg[0] != '-') {
      ReducedImageRef im = ReducedImageNew(arg);
      if (!im || !im->IsValid()) { 
	cerr << args[0] << ": " << arg << " is not a valid dbimage\n";
	continue;
      }
      imList.push_back(im);
      continue;
    }

    switch (arg[1]) {
      //case 'b': { backMethod = (BackgroundMethod) atoi(args[++i]); continue; }
    case 'c': { stackMethod = (StackingMethod) atoi(args[++i]); continue; }
    case 'h':  usage(args[0]);
    case 'w': { weightMethod = (WeightingMethod) atoi(args[++i]); continue; }
    case 'p': { scaleMethod = (PhotoScalingMethod) atoi(args[++i]); continue; }
    case 'o': { outName = args[++i]; continue; }
    case 'z': { phoRef = args[++i]; continue; }
    default: usage(args[0]);
    }
  }

  if (imList.empty()) { cerr << args[0] << ": no images to combine\n";  return EXIT_FAILURE; }

  cout << args[0]<< ": " << imList.size() << " images to combine\n"
       << "   - output image    : " << outName << endl
       << "   - photometric ref : " << phoRef << endl
       << "   - scaling method  : " << name_of_scalingMethod(scaleMethod) << endl
    //       << "   - back method     : " << name_of_backgroundMethod(backMethod) << endl
       << "   - combine method  : " << name_of_stackingMethod(stackMethod) << endl
       << "   - weight method   : " << name_of_weightingMethod(weightMethod) << endl;
  
  try {
    ImageSum stack(outName,
		   imList,
		   phoRef,
		   weightMethod,
		   stackMethod,
		   scaleMethod);
    stack.Execute(DoFits | DoCatalog | DoWeight | DoSatur);
  } catch(PolokaException p) {
    p.PrintMessage(cerr);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
