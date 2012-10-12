#include <string>
#include <iostream>
#include <fstream>

#include "allreducedimage.h"
#include "transformedimage.h"
#include "reducedutils.h"
#include "imagesum.h"
#include "polokaexception.h"

static void usage(const char* ProgName)
{
  cerr << ProgName << " [OPTION...] FILE\n"
       << "Shift and stack subtractions with Poloka\n"
       << "  -c METHOD  : pixel combining method:\n"
       << "               1  Weighted average\n"
       << "               2  Clipped weighted average\n"
       << "               3  Median (default)\n"
       << "               4  Adaptive weighted average\n"
       << "  -o DBIMAGE : output stack DBIMAGE (default: substack)\n"
       << "  -r DBIMAGE : geometric reference (default: first of file)\n";
  exit(EXIT_FAILURE);
}

ReducedImageRef subShift(const string& name,
			 const ReducedImageRef ref,  
			 const double& dx, const double& dy) {

  string shiftname = "Shift_" + ref->Name() + "_" + name;
  {
    ReducedImageRef shifted = new ReducedImage(shiftname);
    if (shifted && shifted->HasImage() && shifted->HasWeight())
      return shifted;
  }
  ReducedImage im(name);
  if (!im.IsValid()) return ReducedImageRef();
  
  GtransfoRef refToIm = FindTransfoFromWCS(*ref, im);
  GtransfoRef shift = new GtransfoLinShift(-dx, -dy);
  GtransfoRef refToImShifted = GtransfoCompose(shift, refToIm);
  GtransfoRef imToRef = FindTransfoFromWCS(im, *ref);
  GtransfoRef shiftinv = new GtransfoLinShift(dx, dy);
  GtransfoRef imShiftedToRef = GtransfoCompose(imToRef, shiftinv);

  ImageGtransfo *imTransfo = new ImageGtransfo(refToImShifted,
					       imShiftedToRef,
					       ref->UsablePart(),
					       ref->Name());

  TransformedImage imResampled(shiftname, im, imTransfo);
  imResampled.Execute(DoFits | DoWeight);
  delete imTransfo;
  return imResampled.Clone();
}

int main(int nargs, char **args) {

  if (nargs < 2) usage(args[0]);

  string outName = "substack";
  const char* fileName;
  StackingMethod combMethod = Median;
  ReducedImageRef ref;

  for (int i=1; i<nargs; ++i) {
    char *arg = args[i];
    if (arg[0] == '-') {
      switch (arg[1]) {
      case 'c': { combMethod = (StackingMethod) atoi(args[++i]); continue; }
      case 'r': { 
	ref = ReducedImageNew(args[++i]);
	if (ref && ref->IsValid())
	  continue;
	else {
	  cerr << args[0] << ": reference " << args[i] << " invalid\n";
	  return EXIT_FAILURE;
	}
      }
      case 'o': outName = args[++i]; break;
      default: usage(args[0]);
      }
      continue;
    }
    fileName = arg;
  }
    
  ReducedImageList subList;

  ifstream ifs(fileName);
  if (ifs.is_open()) {
    char c;
    while (ifs >> c) {
      ifs.unget();
      string name; double mjd, dx, dy;
      ifs >> name >> mjd >> dx >> dy;
      if (ref) {
	ReducedImageRef sub = subShift(name, ref, dx, dy);
	if (sub && sub->IsValid())
	  subList.push_back(sub);
	else
	  cerr << args[0] << ": not a valid dbimage: " << name << endl;
      } else {
	ref = ReducedImageNew(name);
	subList.push_back(ref);
      }
    }
    ifs.close();
  }
  
  if (subList.empty()) {
    cerr << args[0] << ": no valid input images provided\n";
    return EXIT_FAILURE;
  }

  if (outName.empty()) {
    cerr << args[0] << ": no output name provided\n";
    return EXIT_FAILURE;
  }
  
  try {
    ImageSum stack(outName,
		   subList,
		   ref->Name(),
		   PointSourceOptimal,
		   combMethod,
		   ZeroPointDiff);
    stack.Execute(DoFits | DoWeight);
    
  } catch(PolokaException p) {
    p.PrintMessage(cerr);
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
