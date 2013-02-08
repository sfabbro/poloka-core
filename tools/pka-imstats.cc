#include <iostream>
#include <string>

#include <poloka/dbimage.h>
#include <poloka/fitsimage.h>
#include <poloka/frame.h>
#include <poloka/fileutils.h>

static void usage(const char* progname) {
  cerr << "Usage: " << progname << " DBIMAGE...\n"
       << "Print sky, rms for image and weight*image\n";
  exit(EXIT_FAILURE);
}

int main(int argc,char **args)
{
  if (argc <=1) usage(args[0]);

  bool ok = true;
  for (int i=1; i< argc; ++i)
    {
      DbImage dbimage(args[i]);
      if (!dbimage.IsValid())
	{
	  cerr << args[0] << ": invalid file : "  << args[i] << endl;
	  ok = false;
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
      if (im && w)
	{
	  if (FileExists(dbimage.FitsSaturName()))
	    {
	      FitsImage s(dbimage.FitsSaturName());
	      *w *= 1-s;
	    }
	  im_w_stat = ImageAndWeightError(*im, *w);
	  cout << args[i] << ' ' 
	       << " im: ( " << mean_im << ","<< sigma_im << ")"
	       << " w: ("   << mean_w << "," << sigma_w << ")"
	       << " sig(im*sqrt(w)): " << im_w_stat 
	       << endl;
	}
      else
	{
	  cerr << args[0] << ": " << args[i] << " is missing either image or weight\n";
	  ok = false;
	}
      if (im) delete im;
      if (w)  delete w;
    }

  return ok? EXIT_SUCCESS : EXIT_FAILURE;
}
