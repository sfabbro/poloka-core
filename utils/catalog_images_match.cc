#include "reducedimage.h"
#include "dicstar.h"
#include "fitsimage.h"
#include "wcsutils.h"
#include "stringlist.h"
#include "gtransfo.h"


#include "string"

using namespace std;


#include "fstream"
static void read_list(const char * FileName, StringList &Names)
{
  if (!FileExists(FileName))
    {
      cout << " don't find " << FileName << endl;
      exit (-1);
    }
  ifstream f(FileName);
  string name;
  while (f >> name) Names.push_back(name);
  f.close();
}

static void usage(const char *prog)
{
  cout << " usage " << endl;
  cout << prog << " -c <catalog name> dbimages ... " << endl;
  cout << " or :" << endl;
  cout << prog << " -c <catalog name> -l <name of a file that contains dbimage names> " << endl;
  exit (-1);
}
  
#include "imageutils.h"

int main(int nargs, char **args)
{
  StringList imageNames;
  string catalogName;
  for (int i=1; i < nargs; ++i)
    {
      char *arg = args[i];
      if (arg[0] != '-')
	{
	  imageNames.push_back(args[i]);
	  continue;
	}
      else 
	{
	  switch (arg[1])
	    {
	    case 'c' : catalogName = args[++i]; break;
	    case 'l' : read_list(args[++i], imageNames); break;
	    default : 
	      {
		cout << " don't understand " << arg << endl;
		usage(args[0]);
	      }
	    }
	}
    }// end loop on args

  if (catalogName == "")
    {
      cout << " you have to provide a catalog " << endl;
      usage(args[0]);
    }

  DicStarList catalog(catalogName);

  for (StringCIterator it=imageNames.begin(); it != imageNames.end(); ++it)
    {
      ReducedImage im(*it);
      string fitsName = im.FitsName();
      if (!FileExists(fitsName))
	{
	  cout << " Image " << im.Name() 
	       << " has not associated fits header " << endl;
	  continue;
	}
      FitsHeader head(fitsName);
      Gtransfo *wcs;
      if (!WCSFromHeader(head, wcs))
	{
	  cout << " no wcs fo image " << im.Name() << endl;
	  continue;
	}
      Frame imageFrame(head);
      imageFrame.CutMargin(-50); //actually enlarges the frame
      Frame raDecImageFrame = ApplyTransfo(imageFrame, *wcs, LargeFrame);
      DicStarList inImage;
      catalog.ExtractInFrame(inImage,raDecImageFrame);
      cout << " ***** " << im.Name() << " contains " << inImage.size() << " objects " << endl;
      delete wcs;
    }

}

