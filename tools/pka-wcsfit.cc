#include <iostream>

#include <poloka/fitsimage.h>
#include <poloka/dbimage.h>
#include <poloka/dbconfigexception.h>
#include <poloka/fileutils.h>
#include <poloka/usnoutils.h>
#include <poloka/polokaexception.h>

static void usage(const char *progname)
{
  cerr << "Usage: " << progname << "[OPTIONS]... DBIMAGE...\n"
       << "Compute a WCS header for a image\n\n"
       << "    -a CATALOG: use CATALOG as reference to match stars (ra dec mag), superseeds config value\n"
       << "    -c CONFIG : use CONFIG file for options\n"
       << "    -i FITS   : use a FITS image instead of a DBIMAGE (with -l option)\n"
       << "    -l LIST   : use a catalog LIST instead of default from DBIMAGE (with -i option)\n"
       << "    -n        : just print the transfo information, won't alter WCS header\n"
       << "    -o        : overwrite previous WCS\n\n";
  exit(EXIT_FAILURE);
}


int main(int argc, char **argv)
{
  if (argc < 2) usage(argv[0]);
  list<string> dbImageList;

  char *fitsName = NULL;
  char *listName = NULL;
  bool overwrite = false;
  for (int i=1;  i<argc; ++i) // first loop on arguments to read datacards name, if any
    {
      char *arg = argv[i];
      if (arg[0] != '-') continue;
      if (arg[1] == 'c')
	{
	  MatchPrefs.ReadCards(argv[++i]);
	  break;
	}
    }
  for (int i=1;  i<argc; ++i)
    {
      char *arg = argv[i];
      if (arg[0] != '-') 
	{
	  dbImageList.push_back(arg);
	  continue;
	}
      switch (arg[1])
	{
	case 'i' : fitsName = argv[++i]; break;
	case 'l' : listName = argv[++i]; break; 
	case 'n' : MatchPrefs.writeWCS = false; break;
	case 'o' : overwrite = true; break;
	case 'c' : i++; break; // datacards were already read above
	case 'a' : i++;
	  try
	    {
	      MatchPrefs.astromCatalogName =
		DbConfigFindCatalog(argv[i]);
	    }
	  catch (DbConfigException &e)
	    {
	      cerr << e.message() << endl
		   << argv[0] << ": could not locate your astrom catalog\n";
	      return EXIT_FAILURE;
	    }
	  break;
	default: 
	  cerr << argv[0] << ": don't understand argument " << argv[i] << endl;
	  return EXIT_FAILURE;
	}
    }

  if (fitsName && listName)
    {

      try {

	FitsHeader head(fitsName);
	if (head.HasKey("DZEROUSN")  && !overwrite && MatchPrefs.writeWCS && !MatchPrefs.asciiWCS)
	  {
	    cout << argv[0] << ": WCS match already done for " << fitsName << endl;
	    return EXIT_SUCCESS;
	  }	
	if (UsnoProcess(fitsName,listName, NULL))
	  return EXIT_SUCCESS;
	else
	  return EXIT_FAILURE;
      } catch(PolokaException p) {
	p.PrintMessage(cout);
	return EXIT_FAILURE;
      }
    }
  
  
  bool ok = true;
  for (list<string>::iterator i=dbImageList.begin(); i!= dbImageList.end(); ++i)
    {
      string name = *i;
      try {

      DbImage dbimage(name);
      if (!dbimage.IsValid())
	{
	  cerr << argv[0] << ": " << name << " is not a valid DbImage\n";
	  ok = false;
	  continue;
	}

      string catalogName = dbimage.ImageCatalogName(SExtractor);

      string fitsFileName = dbimage.FitsImageName(Calibrated);

      {
	FitsHeader head(fitsFileName);
	if (head.HasKey("DZEROUSN") && !overwrite && MatchPrefs.writeWCS && !MatchPrefs.asciiWCS)
	  {
	    cout << argv[0] << ": WCS match already done for " << name << endl;
	    continue;
	  }	
      }
      ok = UsnoProcess(fitsFileName, catalogName, &dbimage);
      continue;

      } catch(PolokaException p) {
	p.PrintMessage(cout);
	ok = false;
      }
    } /* end of loop on arguments */
  
  return ok? EXIT_SUCCESS : EXIT_FAILURE;
}
