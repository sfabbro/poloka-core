#include <iostream>
#include <string>

#include "dbimage.h"
#include "fitsimage.h"
#include "frame.h"
#include "fileutils.h"

int main(int argc,char **args)
{
if (argc <=1)
  {
    cout << " skylev <dbimage(s)> " << endl;
    exit(-1);
  }
  for (int i=1; i< argc; ++i)
    {
    DbImage dbimage(args[i]);
    if (!dbimage.IsValid())
      {
	cerr << args[0] << " : invalid file : "  << args[i] << endl;
        continue;
      }
    Pixel mean_im = 0, sigma_im = -1;
    FitsImage *im = NULL;
    if (FileExists(dbimage.FitsImageName(Calibrated)))
      {
	im = new FitsImage(dbimage.FitsImageName(Calibrated));
	FitsHeader &head = *im;
	im->SkyLevel(Frame(head),& mean_im, &sigma_im);
      }
    Pixel mean_w = 0, sigma_w = -1;
    FitsImage *w = NULL;
    if (FileExists(dbimage.FitsWeightName()))
      {
	w = new FitsImage(dbimage.FitsWeightName());
	FitsHeader &head = *im;
	w->SkyLevel(Frame(head),& mean_w, &sigma_w);
      }
    double im_w_stat = -1;
    if (im&&w)
      {
	im_w_stat = ImageAndWeightError(*im,*w);
	cout << args[i] << ' ' 
	     << " im: ( " << mean_im << ","<< sigma_im << ")"
	     << " w: ("   << mean_w << "," << sigma_w << ")"
	     << " sig(im*sqrt(w)): " << im_w_stat 
	     << endl;
      }
    else
      cout << " miss either image or weight " << endl;
    if (im) delete im;
    if (w)  delete w;
    }
  return 1;
}
