#include <iostream>
#include <stdlib.h> // for strtod
#include <math.h> // for fabs 
#include <iomanip> // for setprecision

#include "fitsimage.h"
#include "wcsutils.h"
#include "astroutils.h"
#include "fileutils.h"
#include "gtransfo.h"
#include "basestar.h"
#include "frame.h"

void usage()
{
  cout 
    << " astrom <image> -p x y [-d][-v] translates coordinates of x,y in ra,dec" << endl
    << " astrom <image> -p @filename [-d] translates coordinates of x,y in ra,dec" << endl
    << "or" << endl 
    << " astrom <image> -c alpha delta [-d] translates coordinates  alpha delta into pixels" << endl
    << " option -d implies decimal degrees" << endl
    << " option -v dumps WCS" << endl;
  cout << " note : pixel coordinates definition : leftmost bottom pixel center = (1.,1.) " << endl;
  cout << " Hence coordinates from toads catalogs (se.list) should be shifted by 1"<< endl;
}


static void convert_direct(const Gtransfo *Pix2RaDec, 
			   const double x, const double y,
			   const bool decimal,
			   const char* Remainder = NULL)
{
  double ra,dec;
  Pix2RaDec->apply(x,y,ra,dec);
  if (decimal)
    cout << setprecision(10) << ra << " " << dec;
  else
    cout << RaDegToString(ra) << " " << DecDegToString(dec);
  if (Remainder) cout << " " << Remainder;
  else cout << endl;
}

struct Coordinates {
  double x,y;
  string ras, decs;

  // put a default constructor to enventually trace uninitialized data
  Coordinates() : x(-1), y(-1), ras("-1"), decs("92") {};

};


int main(int argc, char **argv)
{
  if (argc < 4)
    {
      usage();
      return 1;
    }
  string fitsName;
  double x;
  double y;
  bool pix2angles = false;
  bool angles2pix = false;
  bool decimal = false;
  bool verbose = false;
  vector<Coordinates> to_convert;
  FILE *list = NULL;
  int status = EXIT_SUCCESS;
  for (int i=1; i< argc; ++i)
    {
      char *arg = argv[i];
      
      if (arg[0] != '-')
	{
	  if (fitsName != "") 
	    {
	      cout << " do not understand " << arg << endl;
	      usage();
	      return EXIT_FAILURE;
	    }
	  fitsName = arg; continue;
	}
      switch (arg[1])
	{
	case 'p' : 
	  pix2angles = true;
	  i++; 
	  if (argv[i][0] =='@')
	    {
	      const char* fileName = argv[i]+1; 
	      list = fopen(fileName,"r");
	      if (!list) 
		{cerr << " cannot open " << fileName << endl; exit(-1);}
	      break;
	    }
	  else 
	    {
	      if (angles2pix)
		{
		  cerr << "cannot mix convertions both ways, choose -p or -c " << std::endl;
		  return EXIT_FAILURE;
		}
	      Coordinates c;
	      c.x = atof(argv[i]) - MEMPIX2DISK; 
	      i++; c.y = atof(argv[i]) - MEMPIX2DISK; 
	      to_convert.push_back(c);
	    }
	  break;
	case 'd' : 
	  decimal = true;
	  break;
	case 'c' :
	  {
	    if (pix2angles)
	      {
		cerr << "cannot mix convertions both ways, choose -p or -c " << std::endl;
		return EXIT_FAILURE;
	      }
	    angles2pix = true;
	    Coordinates c;
	    i++; c.ras  = argv[i];
	    i++; c.decs = argv[i];
	    to_convert.push_back(c);
	    break;
	  }
	case 'v' :
	  verbose = true;
	  break;
	default: 
	  cerr << " astrom : bad arg " << argv[i] << endl;
	  return EXIT_FAILURE;
	}
    }// end looping on args

  if (!FileExists(fitsName))
    {
      cerr << " cannot open " << fitsName << endl;
      return 0;
    }
  Gtransfo *Pix2RaDec;
  FitsHeader head(fitsName);
  if (!WCSFromHeader(head, Pix2RaDec))
    {
      cerr << " do not find the expected WCS in " << fitsName  << endl;
      return 0;
    }
  if (verbose) cout << " WCS : " << endl << *Pix2RaDec << " ======= " << endl; 


  if (pix2angles)
    {
      if (!list) 
	{
	  for (unsigned k=0; k < to_convert.size(); ++k)
	    {
	      Coordinates &c = to_convert[k]; 
	      convert_direct(Pix2RaDec, c.x, c.y, decimal);
	    }
	}
      else
	{
	  char line[2048];
	  char* remainder = NULL;
	  while(fgets(line,2048,list))
	    {
	      if (line[0] == '#' || line[0] == '@') {cout << line; continue;}
	      char *next_stuff = line;
	      x = strtod(line, &next_stuff);
	      if (next_stuff != line) y = strtod(next_stuff, &remainder);
	      if (next_stuff == remainder)
		{
		  cerr << " ERROR : cannot decode 2 coordinates in :" << line;
		  status = EXIT_FAILURE;
		  break;
		}
	      convert_direct(Pix2RaDec,x,y,decimal, remainder);
	    }
	  fclose(list);
	}
    }
  else  if (angles2pix)
    {
      // 0.01 is the precision ot invertion (in pixels);
      Gtransfo *RaDec2Pix = Pix2RaDec->InverseTransfo(0.01, Frame(head));
      for (size_t i=0; i<to_convert.size(); ++i)
	{
	  Coordinates &c = to_convert[i];
	  double ra = RaStringToDeg(c.ras); // in principe decodes both sexagesimal and decimal.
	  double dec = DecStringToDeg(c.decs); // same comment.
	  double x,y;
	  RaDec2Pix->apply(ra,dec,x,y);
	  cout << x +MEMPIX2DISK << " " << y + MEMPIX2DISK << endl;
	  // check that the result makes sense:
	  double raf, decf;
	  Pix2RaDec->apply(x,y,raf,decf);
	  double error = fabs(ra-raf)/cos(M_PI*dec/180.)+fabs(dec-decf);
	  if (error > 0.2/3600.)
	    {
	      cerr << " large error when inverting WCS " << endl;
	      status = EXIT_FAILURE;
	    }
	}
    }
  return status;
}
