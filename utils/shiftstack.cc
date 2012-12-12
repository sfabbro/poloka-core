#include <string>
#include <iostream>
#include <fstream>

#include "allreducedimage.h"
#include "transformedimage.h"
#include "basestar.h"
#include "daophotio.h"
#include "reducedutils.h"
#include "imagesum.h"
#include "polokaexception.h"

static void usage(const char* ProgName)
{
  cerr << ProgName << " [OPTION...] FILE\n"
       << "Shift and stack DbImage\n"
       << "FILE line format: DbImage ra_shift dec_shift\n"
       << "OPTIONS:\n"
       << "  -c METHOD  : pixel combining method:\n"
       << "               1  Weighted average\n"
       << "               2  Clipped weighted average\n"
       << "               3  Median (default)\n"
       << "               4  Adaptive weighted average\n"
       << "  -o DBIMAGE : output stack DBIMAGE (default: shiftstack)\n"
       << "  -r DBIMAGE : geometric reference (default: first of file)\n"
       << "  -l STRING  : use star list STRING to compute transfos instead of shifts\n";
  exit(EXIT_FAILURE);
}

static BaseStarList imStarList(const string& slist, const ReducedImage& im) {
  string suf = FileExtension(slist);
  BaseStarList bList;
  if (suf == "list") { // poloka shit
    bList.read(im.Dir() + slist);
  } else { // daophot shit
    DaoStarList imStars;
    ReadDaoList(im.Dir() + slist, imStars);
    if (imStars.empty()) {
      cerr <<  im.Name()  << " : empty " << slist << endl;
      return bList;
    }
    bList = Dao2Base(imStars);
  }
  return bList;
}

ReducedImageRef imShift(const string& name,
			const ReducedImageRef ref,  
			const double& dra,
			const double& ddec,
			const string& slist) {

  string shiftname = "ShiftStack_" + ref->Name() + "_" + name;
  {
    ReducedImageRef shifted = ReducedImageNew(shiftname);
    if (shifted && shifted->HasImage() && shifted->HasWeight())
      return shifted;
  }

  ReducedImage im(name);
  if (!im.IsValid()) return ReducedImageRef();

  GtransfoRef imToRefShift, shiftRefToIm;

  if (slist.empty()) {
    GtransfoRef imToRef = FindTransfo(im, *ref);
    GtransfoRef refToIm = FindTransfo(*ref, im);
    GtransfoRef raDecToPix = ref->RaDecToPixels();
    double rac = ref->RaDeg2000();
    double dec = ref->DecDeg2000();
    double x, y, dx, dy;
    raDecToPix->apply(rac, dec, x, y);
    raDecToPix->apply(rac + dra/cos(M_PI*dec/180.), dec + ddec, dx, dy);
    dx -= x;
    dy -= y;
    GtransfoRef shift = new GtransfoLinShift(dx, dy);
    GtransfoRef shiftinv = new GtransfoLinShift(-dx, -dy);
    imToRefShift = GtransfoCompose(shift, imToRef);
    shiftRefToIm = GtransfoCompose(refToIm, shiftinv);
  } else {
    BaseStarList refList = imStarList(slist, *ref);
    BaseStarList imList = imStarList(slist, im);
    imToRefShift = ListMatch(imList, refList);
    shiftRefToIm = ListMatch(refList, imList);
  } 

  ImageGtransfo *imTransfo = new ImageGtransfo(shiftRefToIm,
					       imToRefShift,
					       ref->UsablePart(),
					       ref->Name());

  TransformedImage imResampled(shiftname, im, imTransfo);
  imResampled.Execute(DoFits | DoWeight);
  delete imTransfo;
  return imResampled.Clone();
}

int main(int nargs, char **args) {

  if (nargs < 2) usage(args[0]);

  string outName = "shiftstack";
  string listName;
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
      case 'o': outName  = args[++i]; break;
      case 'l': listName = args[++i]; break;
      default: usage(args[0]);
      }
      continue;
    }
    fileName = arg;
  }
    
  ReducedImageList imList;

  ifstream ifs(fileName);
  if (ifs.is_open()) {
    char c;
    while (ifs >> c) {
      ifs.unget();
      string line;
      getline(ifs, line);
      if (line.empty() || line[0] == '#') continue;
      istringstream iline(line);
      string name; double dra, ddec;
      iline >> name >> dra >> ddec;
      if (ref) {
	ReducedImageRef im = imShift(name, ref, dra, ddec, listName);
	if (im && im->IsValid())
	  imList.push_back(im);
	else
	  cerr << args[0] << ": not a valid dbimage: " << name << endl;
      } else {
	ref = ReducedImageNew(name);
	imList.push_back(ref);
      }
    }
    ifs.close();
  }
  
  if (imList.empty()) {
    cerr << args[0] << ": no valid input images provided\n";
    return EXIT_FAILURE;
  }

  if (outName.empty()) {
    cerr << args[0] << ": no output name provided\n";
    return EXIT_FAILURE;
  }
  
  try {
    ImageSum stack(outName,
		   imList,
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
