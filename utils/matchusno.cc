#include <iostream>

#include "fitsimage.h"
#include "dbimage.h"
#include "fileutils.h"
#include "usnoutils.h"

static void usage(const char *pgname)
{
  cout << pgname << ' ' << " <dbimage name .... > (for matching with usno)" << endl;
  cout <<   " or " << pgname << ' ' << " -i <fitsImage> -l <cat> [-n](just print)" << endl;
  cout << " [-o] to overwrite previous match " << endl;
}


int main(int argc, char **argv)
{
  list<string> dbImageList;
  if (argc < 2)
    {
      usage(argv[0]);  return 0;
    }

  char *fitsName = NULL;
  char *listName = NULL;
  bool write = true;
  bool overwrite = false;
  for (int i=1;  i<argc; ++i)
    {
      char *arg = argv[i];
      if (arg[0] != '-') 
	{
	  dbImageList.push_back(string(arg));
	  continue;
	}
      switch (arg[1]) 
	{
	case 'i' : i++; fitsName = argv[i]; break;
	case 'l' : i++; listName = argv[i]; break; 
	case 'n' : write = false; break;
	case 'o' : overwrite = true; break;
	default : usage(argv[0]); exit(0);
	}
    }

  if (fitsName && listName)
    {
      FitsHeader head(fitsName);
      if (head.HasKey("DZEROUSN")  && !overwrite && write)
	{
	  cout << "Match already done! " << endl;
	  exit(1);
	}	
      UsnoProcess(fitsName,listName, NULL, write);
      return 1;
    }
  

  for (list<string>::iterator i=dbImageList.begin(); i!= dbImageList.end(); ++i)
    {
      string name = *i;
      DbImage dbimage(name);
      if (!dbimage.IsValid())
	{
	  cerr << " Be careful ! " << name << " must be an image name ! " << endl;
	  continue;
	}

      string catalogName = dbimage.ImageCatalogName(SExtractor);
      if (!FileExists(catalogName.c_str()))
	{
	  cerr << "The SExtractor Catalogue associated to image " << name << " doesn't exist !! " << endl;
	  continue;
	}  
      string fitsFileName = dbimage.FitsImageName(Calibrated);
      if (!FileExists(fitsFileName))
	{
	  cerr << "The  calibrated fits image " << name << " doesn't exist !! " << endl;
	  continue;
	}
      {
	FitsHeader head(fitsFileName);
	if (head.HasKey("DZEROUSN") && !overwrite && write)
	  {
	    cout << "Match already done for " << name << endl;
	    continue;
	  }	
      }
      UsnoProcess(fitsFileName, catalogName, &dbimage,write);

#ifdef STORAGE
      cout << "The file \"" << MatchFile << "\" is filled" << endl; 
      char *fluxfile = FluxFile(Flux_x, Flux_y, &dbimage); 
      cout << "The file \"" << fluxfile << "\" is filled" << endl; 

      // The following allow to plot the straight line fit on the data 
      double x_min = -2.5*log10(Flux_x.Max());
      double x_max = -2.5*log10(Flux_x.Min());
      //cout << x_min << " " << x_max << " " << Flux_x.Size() << endl;
      char *DroiteFile = StraightLineFile(1., zeropoint, x_min, x_max, &dbimage);
   
      cout << "The File \"" << DroiteFile << "\" is filled" << endl;
      char *kumacfile = DroiteKumac(&dbimage, MatchFile, fluxfile, DroiteFile);
      cout << "The Kumac File \"" << kumacfile << "\" is filled" << endl;
      cout << endl << " ****************************** " << endl << endl; 
#endif
    } /* end of loop on arguments */

  return 1;
}


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

