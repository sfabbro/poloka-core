#include "fitsimage.h"
#include "frame.h"
#include "fitstoad.h"
#include "polokaexception.h"

static void usage()
{
  cout << " trim -t '[x0:x1,y0:y1]' <fitsfile(s)> \n" 
       << "      will cut the image to save the illuminated frame \n"
       << "      -t '[x0:x1,y0:y1]' : define the section to cut in FITS coordinates\n"
       << " 	             (if omitted will take the one in header)\n"
       << endl;
  
  exit(-1);
}

int main(int argc,char **args)
{

  if (argc <=1) usage();

  // keep it fast, so avoid list<string>
  bool use_header = true;
  int istart = 1;
  Frame illu;
  
  char *arg = args[1];
  if (arg[0] == '-')
    {
      if (arg[1] != 't') usage();
      int x0,x1,y0,y1;
      if (sscanf(args[2],"[%d:%d,%d:%d]",&x0,&x1,&y0,&y1) == 4)
	{
	  if (x0>x1) swap(x0,x1);
	  if (y0>y1) swap(y0,y1);
	  illu = Frame(Point(x0-1,y0-1), Point(x1-1, y1-1));
	  use_header = false;
	}
      else cerr << " cannot decode 4 integers in " << args[2] << endl;
      istart += 2;
    }

  bool ok = true;
  for (int i=istart; i< argc; ++i)
    {
      try{
      FitsImage image(args[i],RW);
      image.EnableWrite(false); // hold it
      if (!image.IsValid())
	{
	  cerr << " trim : invalid fits file : "  << args[i] << endl;
	  continue;
	}
      if (use_header) illu = TotalIlluRegion(image);
      if (image.Trim(illu)) image.EnableWrite(true);
      // keep the original scaling : 
      image.Write( /*force_bscale = */ true);
      }catch(PolokaException p){
	p.PrintMessage(cout);
	ok = false;
      }
    }
  if(ok)
   return EXIT_SUCCESS;
  return EXIT_FAILURE;
}


