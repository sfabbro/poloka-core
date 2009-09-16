#include "lightcurvefile.h"
#include "polokaexception.h"

#include <iostream>
#include <string>

using namespace std;

#define DEF_LCFILE "makelightcurve.conf"

static void usage(const char *prog)
{
  std::cout << "usage : " << endl;
  std::cout << prog << "[<lightcurve_file=\""DEF_LCFILE"\">]" 
	    << " [-v](write vignettes) " << endl
            << " [-m}(write matrices) "  << endl
	    << " [-d](one output directory per object) " << endl
	    << " [-c] <calibrationCatalog>(fits objects in catalog on images in lightcurve_file) " << endl;
  exit(-1);
}


int main(int nargs, char **args)
{

  string filename(DEF_LCFILE);
  string calibrationCatalog;
  bool writeVignettes = false;
  bool writeMatrices = false;
  bool oneDirPerObj = false;
  for (int i=1; i< nargs; ++i)
    {
      char *arg = args[i];
      if (arg[0] == '-')
	switch(arg[1])
	  {
	  case 'h' : usage(args[0]); break;
	  case 'v' : writeVignettes = true; break;
	  case 'm' : writeMatrices = true; break;  
	  case 'd' : oneDirPerObj = true; break;
	  case 'c' : calibrationCatalog=args[++i];continue;break;
	  default : std::cout << " don't understand " << arg << endl;
	    usage(args[0]);
	  }
      else
	filename = arg;
    }
  bool success = true;
  try 
    {
      LightCurveFile lcf(filename);
      if (writeVignettes) lcf.PleaseWriteVignettes();
      if (writeMatrices) lcf.PleaseWriteMatrices();
      if (oneDirPerObj) lcf.PleaseOneDirPerObject();
      if (calibrationCatalog != "")
	{
	  success = lcf.SimPhotFitAllCalib(calibrationCatalog);
	}
      else success = lcf.SimPhotFitAll();
    }
  catch (PolokaException e)
    {
      cout << e.message() << endl;
      success = false;
    }
  return ((success)? EXIT_SUCCESS :  EXIT_FAILURE) ;
}
  
      
