#include <iostream>
#include <iomanip>
#include "reducedimage.h"
#include "fileutils.h"
#include "fitsimage.h"

void usage(char *progName)
{
  cerr << progName << " [-s sort by seeing ] [-d sort by date (default)] <DbImages> " << endl;
}

int main(int nargs, char **args)
{
  if (nargs <= 1){usage(args[0]); exit(1);}
  cout << endl << " Checking images " << endl;
  int old = cout.precision();
  long f = cout.flags() ;
  cout << resetiosflags(ios::scientific);
  cout << setiosflags(ios::fixed);
  ReducedImageList alls;
  unsigned maxs = 0;
  bool seeingsort = false;
  bool datesort = true;
  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      if (arg[0] != '-') 
	{
	  string name = string(args[i]);       	    
	  ReducedImage *current = new ReducedImage(name);
	  if (FileExists(current->FitsImageName(Calibrated))) 
	    {
	      if (name.length() > maxs) maxs = name.length();
	      alls.push_back(current);
	    }
	  else delete current;
	  continue;
	}
      switch (arg[1])
	{
	case 's' : seeingsort = true; datesort = false; break;
	case 'd' : seeingsort = false; datesort = true; break;
	}

    }

  if (seeingsort) sort(alls.begin(), alls.end(), IncreasingSeeing);
  if (datesort) sort(alls.begin(), alls.end(), IncreasingJulianDate);

  cout << "---------------------------------------------------------------------------------------------------------" << endl;
  cout << setw(maxs+2) << setiosflags(ios::left)<< "Image" << setiosflags(ios::right);
  cout << setw(10) << "Date" << setw(5) << "Band" << setw(7) << "Exp." << setw(7)
       << "FWHM" << setw(10) << "Sky" << setw(9) << "Satur" << setw(5) 
       << "ZP" << setw(12) << "Target" << endl;
  cout << "---------------------------------------------------------------------------------------------------------" << endl;
  for (unsigned int i=0; i<alls.size(); ++i)
    {
      ReducedImage *current = alls[i];
      cout << setiosflags(ios::left);      
      cout << setw(maxs+2) << current->Name().c_str();
      cout << setiosflags(ios::right);      
      cout << setw(12) << current->Date().c_str() << setw(3) << current->Band().c_str();
      cout << setw(7) << setprecision(0) << current->Exposure();
     
      if (!current->ActuallyReduced()) 
	{
	  cout << "  *******NOT REDUCED******* " << endl; 
	  continue;
	}
      cout << setw(7) << setprecision(2) << current->Seeing()*2.3548*current->PixelSize()
	   << setw(8) << setprecision(0) << current->BackLevel()
	   << setw(8) << setprecision(0) << current->Saturation()
	   << setw(5) << setprecision(1) << current->ZeroPoint();
      cout << "  " << current->Target() << endl;
    }
  cout << setprecision(old); 
  //cout.flags(f);
  return EXIT_SUCCESS;
}
