#include "lightcurvefile.h"
#include "polokaexception.h"

#include <iostream>
#include <string>

using namespace std;

#define DEF_LCFILE "mklc_wnr.conf"

static void usage(const char *prog)
{
  std::cout << "usage : " << endl;
  std::cout << prog << "[<lightcurve_file=\""DEF_LCFILE"\">]" 
	    << " [-v](write vignettes) " << endl
            << " [-m}(write matrices) "  << endl
            << " [-N] n max star to be fitted "  << endl
            << " [-f] fittype : 0=gal+flux 1=flux (default is 0)"  << endl
	    << " [-d](one output directory per object) " << endl 
	    << " [-o] <outputcatalog>  " << endl
	    << " [-c] <calibrationCatalog> " << endl;

  exit(-1);
}


int main(int nargs, char **args)
{

  string filename(DEF_LCFILE);
  string outputCatalog;
  string calibrationCatalog;
  bool writeVignettes = false;
  bool writeMatrices = false;
  bool oneDirPerObj = false;
  int Nmax = -1 ;
  int fit_type=0;
  for (int i=1; i< nargs; ++i)
    {
      char *arg = args[i];
      if (arg[0] == '-')
	switch(arg[1])
	  {
	  case 'h' : usage(args[0]); break;
	  case 'v' : writeVignettes = true; break;
	  case 'N' : Nmax = atoi(args[++i]); break;  
	  case 'm' : writeMatrices = true; break;  
	  case 'd' : oneDirPerObj = true; break;
	  case 'o' : outputCatalog=args[++i];continue;break;
	  case 'c' : calibrationCatalog=args[++i];continue;break;
	  case 'f' : fit_type=atoi(args[++i]);continue;break;
	  default : std::cout << " don't understand " << arg << endl;
	    usage(args[0]);
	  }
      else
	filename = arg;
    }
  cerr << "Fit Type : " << fit_type << endl ;
  bool success = true;
  try 
    {
      LightCurveFile lcf(filename);
      if (writeVignettes) lcf.PleaseWriteVignettes();
      if (writeMatrices) lcf.PleaseWriteMatrices();
      if (oneDirPerObj) lcf.PleaseOneDirPerObject();
      if (calibrationCatalog != "")
	{
	  if (fit_type==0)
	    {
	      cout << " assuming you mean fit_type = 1" << endl;
	      fit_type = 1;
	    }
	  success = lcf.SimPhotFitAllCalib(calibrationCatalog, outputCatalog, fit_type, Nmax);
	}
      else success = lcf.SimPhotFitAll();
    }
  catch (PolokaException e)
    {
      cout << "ERROR : " <<  e.message() << endl;
      success = false;
      abort();
    }
  return ((success)? EXIT_SUCCESS :  EXIT_FAILURE) ;
}
  
      
