#include <subimage.h>
#include <fitsimage.h>
#include <reducedutils.h>
#include <wcsutils.h>

int main(int nargs, char** args) {

  ReducedImage toExtract(args[1]);

  FitsHeader largeHead(toExtract.FitsName());
  Frame largeFrame(largeHead, WholeSizeFrame);
  GtransfoRef largePix2RaDec = WCSFromHeader(largeHead);
  if (!largePix2RaDec)  {
    cerr << "extract_image: cannot handle a large reference without a WCS\n";
    return EXIT_FAILURE;
  }

  GtransfoRef largeRaDec2Pix = 
    largePix2RaDec->InverseTransfo(0.1, largeFrame);
  Frame frame;

  for (int i=2; i<nargs; ++i) {
    // first compute the union of Frame(largeImage)*Frame(anyOtherImage)
    ReducedImage current(args[i]);
    if (current == toExtract) continue;
    FitsHeader currentHead(current.FitsName());
    GtransfoRef current2Large = FindTransfoFromWCS(current, toExtract);
    cout << " extract_image: transfo between " << toExtract.Name() << ' ' << " and " 
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
    cout << " extract_image: after adding " << current.Name() 
	 << ", extract frame = " << frame << endl;
  }
  
  SubImage *subImage = new SubImage("Ext_" + toExtract.Name(), toExtract.Name(), frame);
  subImage->MakeFits();
  subImage->MakeWeight();
  subImage->MakeCatalog();

  return EXIT_SUCCESS;
}
