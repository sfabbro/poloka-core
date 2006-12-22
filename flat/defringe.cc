
#include <iostream>
#include <vector>

#include "image.h"
#include "fitsimage.h"
#include "fitsimagearray.h"
#include "fringeutils.h"
//#include "dbimage.h"
// get cvs version of the code
//#define CVSVERSION "$Revision: 1.2 $"

void DumpHelp(const char *progName) {
  cout     << progName << " removes fringes of a FITS image" << endl
	   << "   usage : " << progName << " <image.fits> -f <fringe> [options]"<< endl
    //<< "       or  " << progName << " <dbimage>    [options]"<< endl
	   << "    [options]:" << endl
	   << "        -nvec #     : number of fringe patterns to use (default is the content of fringe.fits)" << endl
	   << "        -o <file>   : output defringed image (default is overwrite)" << endl
    	   << "        -bg         : substract background" << endl
	   << "        -nsig #     : number of sigma cutoff to compute scalar product (default is 3)" << endl
	   << "        -v          : verbose" << endl
	   << "        --help (-h) : this help" << endl
    	   << endl;
  exit(1);
}


int main(int argc, char**argv) {
  
  if (argc < 2) 
    DumpHelp(argv[0]);
  
  string in_name = ""; 
  string out_name = "";
  string fringefilename = "";
  
  int nvec = 0;
  bool substractbg = false;
  float nsig = 3; // Cuts on number of sigma for the science image
  bool verbose = false;


  for (int i=1; i< argc; i++) {
    if (!strcmp(argv[i],"--help") || !strcmp(argv[i],"-h")) {
      DumpHelp(argv[0]);
    }
    else if (!strcmp(argv[i],"-f")) {
      fringefilename = argv[++i];
      continue;
    }
    else if (!strcmp(argv[i],"-nsig")) {
      nsig = atof(argv[++i]);
      continue;
    }
    else if (!strcmp(argv[i],"-bg")) {
      substractbg = true;
      continue;
    }
    else if (!strcmp(argv[i],"-nvec")) {
      nvec = atoi(argv[++i]);
      continue;
    }
    else if (!strcmp(argv[i],"-o")) {
      out_name = argv[++i];
      continue;
    } 
    else if (!strcmp(argv[i],"-n")) {
      cout << "WARNING: -n is an option of defringe unused here" << endl;
      continue;
    }
    else if (!strcmp(argv[i],"-v")) {
      verbose = true;
      continue;
    }
    in_name = argv[i];
  }


  FitsFileMode fmode = RW;
  if (out_name != "") fmode = RO;
  
  /*
  if (fringefilename == "")  {
    DbImage dbim(in_name);
    in_name = dbim.FitsImageName(Calibrated);
    fringefilename = dbim.FitsFringeName();
  }
  */
  if (fringefilename == "")  {
    cout << "Error with input parameters" << endl;
    cout << "(Try " << argv[0] << " --help)" << endl;
    exit(-1);
  }

  
  if(fmode == RW) {
    FitsImage image(in_name, RW);
    if(!image.IsValid()) {
      cout << "Input image is not valid" << endl;
      cout << "(Try " << argv[0] << " --help)" << endl;
      exit(-1);
    }
    if(verbose)
      cout << "Cut at " << nsig << " sigmas" << endl;
  
    if( FringeUtils::RemoveFringes(image,fringefilename,nvec,nsig,substractbg,verbose))
      return -1;
  }else{
     FitsImage input(in_name, RO);
     FitsImage image(out_name,input,input);
      if(!image.IsValid()) {
      cout << "Input image is not valid" << endl;
      cout << "(Try " << argv[0] << " --help)" << endl;
      exit(-1);
    }
    if(verbose)
      cout << "Cut at " << nsig << " sigmas" << endl;
  
    if( FringeUtils::RemoveFringes(image,fringefilename,nvec,nsig,substractbg,verbose))
      return -1;    
  }

  return 1;
}
