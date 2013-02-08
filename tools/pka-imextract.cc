#include <poloka/subimage.h>
#include <poloka/fitsimage.h>
#include <poloka/reducedutils.h>
#include <poloka/wcsutils.h>

static void usage(const char *progname) {
  cerr << "Usage: " << progname << " REFIMAGE [OPTION] DBIMAGE...\n"
       << "Extract the minimal sub image from REFIMAGE that matches all DBIMAGE\n"
       << "   -o DBIMAGE: output name of the extracted dbimage (default is 'Ext_REFIMAGE')\n";
  exit(EXIT_FAILURE);
}


int main(int nargs, char** args) {

  if (nargs < 2) usage(args[0]);
  ReducedImage toExtract(args[1]);
  string extname = "Ext_" + toExtract.Name();

  FitsHeader largeHead(toExtract.FitsName());
  Frame largeFrame(largeHead, WholeSizeFrame);
  GtransfoRef largePix2RaDec = WCSFromHeader(largeHead);
  if (!largePix2RaDec)  {
    cerr << args[0] << ": cannot handle a large reference without a WCS\n";
    return EXIT_FAILURE;
  }

  GtransfoRef largeRaDec2Pix = 
    largePix2RaDec->InverseTransfo(0.1, largeFrame);
  Frame frame;

  for (int i=2; i<nargs; ++i) {
    char *arg = args[i];
    if (arg[0] == '-') {
      switch (arg[1]) {
      case 'o': { extname = args[++i]; continue; }
      default: usage(args[0]);
      }
    }

    // first compute the union of Frame(largeImage)*Frame(anyOtherImage)
    ReducedImage current(arg);
    if (current == toExtract) continue;
    FitsHeader currentHead(current.FitsName());
    GtransfoRef current2Large = FindTransfoFromWCS(current, toExtract);
    cout << args[0] << ": transfo between " << toExtract.Name() << ' ' << " and " 
	 << current.Name() << endl
	 << *current2Large;

    Frame currentFrame(currentHead, WholeSizeFrame);
    currentFrame = ApplyTransfo(currentFrame, *current2Large);
    /* if we are at the first image in the loop we have to initialize
       the frame, else we just "increment" it so that we have the union
       of all input images at the end */
    if (frame.Area() == 0.)
      frame = (largeFrame*currentFrame);
    else
      frame += (largeFrame*currentFrame);
    cout << args[0] << ": after adding " << current.Name() 
	 << ", extract frame = " << frame << endl;
  }
  // stupid hack to enlarge images to account for distorsions
  frame.CutMargin(-100, -100);

  SubImage subImage(extname, toExtract.Name(), frame);

  bool ok = true;

  if (toExtract.HasImage())
    ok = subImage.MakeFits();
  if (toExtract.HasWeight())
    ok = subImage.MakeWeight();
  if (toExtract.HasSatur())
    ok = subImage.MakeSatur();

  ok = subImage.MakeCatalog();
      
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
