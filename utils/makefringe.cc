#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string>

#include "fileutils.h"
#include "fitsimage.h"
#include "superflat.h"

void usage()
{
  cout << endl <<"makefringe builds a fringe map" << endl
       << " usage : makefringe <fringed flat> <blank flat> [options]" << endl
       << "   <fringed flat>  is usually a superflat FITS file, full of fringes. "<< endl
       << "   <blank flat> is a  master dome or twillight flat FITS file, fringeless." << endl
       << "   [options] are :" << endl
       << "       -f : filter the fringe map with a low pass band filter" << endl
       << "       -o <out name>: Write the fringe map as a FITS file in <out name>"<< endl
       << "       -b : subtract a background image of the constructed fringe map"<< endl
       << "       -c <nsigmas>: Force pixels to be at max or min of <nsigmas>*sigma"<< endl
       << endl;
  
  exit(0);
}


int main(int argc, char**argv)
{
  if (argc < 3) usage();

  bool blankflag = false;
  bool fringedflag = false;
  bool filterflag = false;
  bool writeflag = false;
  bool cutflag = false;
  bool backflag = false;
  string blank_name, fringed_name,fringemap_name;
  double nsigcut = 3.5;
  for (int i=1; i< argc; i++)
    {
      char *arg = argv[i];
      if ((arg[0] != '-') && !(fringedflag)){fringed_name=argv[i];fringedflag=true;}
      if ((arg[0] != '-') && !(blankflag) && (fringedflag)){blank_name=argv[i+1];blankflag=true;}
      if (arg[0] == '-')
	  {
	    switch (arg[1])
	      {
	      case 'f' : filterflag=true;;break;
	      case 'o' : i++;writeflag=true;fringemap_name=argv[i];break;
	      case 'b' : backflag=true;break;
	      case 'c' : i++;cutflag=true;nsigcut=atof(argv[i]);break;
	      default : argc = 3;
	      }
	  }
    }     
  if ( (!blankflag) || (!fringedflag) ) usage();
  
  FitsImage blankFits(blank_name);
  FitsImage fringedFits(fringed_name);

  string blankfilt = blankFits.KeyVal("TOADBAND"); 
  string fringedfilt = fringedFits.KeyVal("TOADBAND"); 
  string blankchip = blankFits.KeyVal("TOADCHIP"); 
  string fringedchip = fringedFits.KeyVal("TOADCHIP"); 

  if (blankfilt != fringedfilt)
    {
      cout << " blank and fringed don't have the same filter ! Stop here !"<<endl;
      return -1;
    }
  if (blankchip != fringedchip)
    {
      cout << " blank and fringed don't have the same chip number ! Stop here !"<<endl;
      return -1;
    }

  
  if (!writeflag) 
    {
      string sccd = blankFits.KeyVal("TOADCHIP");
      fringemap_name = "FringeMap_ccd" + sccd +".fits";
    }
  
  cout << "Making fringe pattern...." << endl;
  Image *fringemapImage = MakeFringePattern(fringedFits, blankFits,nsigcut,backflag,filterflag);
  cout << "done" << endl;


  //FitsImage fringeFits(fringemap_name,blankFits,*fringeImage);
  FitsImage fringeFits(fringemap_name, fringedFits, *fringemapImage);
  fringeFits.AddOrModKey("NUMERATO",fringed_name.c_str(),"Name of FITS file of the used fringed superflat on numerator");
  fringeFits.AddOrModKey("DENOMINA",blank_name.c_str(),"Name of FITS file of the used blank (twillight or dome) flat on denominator");
  //fringeFits.Write(true);
 
   delete fringemapImage;

}
