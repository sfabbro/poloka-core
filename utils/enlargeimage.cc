#include <iostream>
#include <vector>

#include "fitsimage.h"
#include "frame.h"
#include "transformedimage.h"

using namespace std;

static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " dbimage2enlarge newdbimage marginx marginy  " << endl;
  cerr << "  add margins of size marginx (left and right) and marginy (top and bottom)" << endl;
  cerr << "  to image dbimage2enlarge and save it in newdbimage"  << endl;
  exit(EXIT_FAILURE);
}

int main(int nargs, char **args)
{
  // if nothing is given
  if (nargs < 5){usage(args[0]);}
  
  
  CountedRef<ReducedImage> refimage = new ReducedImage(args[1]);
  string newimagename = args[2];
  
  int marginx = atoi(args[3]);
  int marginy = atoi(args[4]);
  
  // fitsimage
  FitsImage reffits(refimage->FitsName());
  
  // create a gtransfo which is just a shift
  CountedRef<Gtransfo>  TransfoFromRef = new GtransfoLinShift(-marginx,-marginy);
  CountedRef<Gtransfo>  TransfoToRef = new GtransfoLinShift(marginx,marginy);
  
  Frame InputImageSize((const Image&)reffits);
  Frame OutputImageSize(InputImageSize.xMin,InputImageSize.yMin,InputImageSize.xMax+2*marginx,InputImageSize.yMax+2*marginy);
  ImageGtransfo transfo(TransfoFromRef,TransfoToRef,OutputImageSize,"");
  
  // now transform image
  TransformedImage transformed(newimagename,*refimage,&transfo);
  transformed.Execute(DoFits | DoCatalog | DoSatur | DoWeight | DoCosmic);

  FitsHeader head(transformed.FitsName());
  head.AddOrModKey("XMARGIN", marginx, "X margin width of enlarged image");
  head.AddOrModKey("YMARGIN", marginy, "Y margin width of enlarged image");
  
  return EXIT_SUCCESS;
}
