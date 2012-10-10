#include <iostream>
#include <sstream>
#include <cstdlib> // for strtod
#include <cmath> // for fabs 
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
    << " astrom <image> -c alpha delta [-s] translates coordinates  alpha delta into pixels" << endl
    << " astrom <image> -c @filename [-s] translates coordinates  alpha delta into pixels" << endl
    << " option -d implies decimal degrees" << endl
    << " option -v dumps WCS" << endl
    << " option -s skips out of frame coordinates" << endl;
  cout << " note : pixel coordinates definition : leftmost bottom pixel center = (1.,1.) " << endl;
  cout << " Hence coordinates from toads catalogs (se.list) should be shifted by 1"<< endl;
}


static void convert_direct(const GtransfoRef Pix2RaDec, 
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
  bool skip_out_of_bounds = false;
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
	case 's' : 
	  skip_out_of_bounds = true;
	  break;
	case 'c' :
	  {
	    angles2pix = true;
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
		if (pix2angles)
		  {
		    cerr << "cannot mix convertions both ways, choose -p or -c " << std::endl;
		    return EXIT_FAILURE;
		  }
	    Coordinates c;
	    c.ras  = argv[i++];
	    c.decs = argv[i];
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
  FitsHeader head(fitsName);
  GtransfoRef Pix2RaDec = WCSFromHeader(head);
  if (!Pix2RaDec)
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
	  char line[8192];
	  char* remainder = NULL;
	  while(fgets(line,8192,list))
	    {
	      if (line[0] == '#')
		{
		  if (strstr(line,"format"))
		    {
		      int format = DecodeFormat(line,"BaseStar");
		      if (format !=2)
			cerr << " WARNING : this code does not convert position errors ... yet. " << endl;
		    }
		  cout << line; continue;
		}
	      if (line[0] == '@') {cout << line; continue;}
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
      Frame frame(head);
      GtransfoRef RaDec2Pix = Pix2RaDec->InverseTransfo(0.01, frame);
      if (!list) 
	{
	  for (size_t i=0; i<to_convert.size(); ++i)
	    {
	      Coordinates &c = to_convert[i];
	      double ra = RaStringToDeg(c.ras); // in principe decodes both sexagesimal and decimal.
	      double dec = DecStringToDeg(c.decs); // same comment.
	      double x,y;
	      RaDec2Pix->apply(ra,dec,x,y);
	      if (skip_out_of_bounds && !frame.InFrame(x,y)) continue;
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
      else
	{
	  char line[8192];
	  string ras,decs;
	  while(fgets(line,8192,list))
	    {
	      if (line[0] == '#' || line[0] == '@') continue;
	      istringstream iline(line);	      
	      if (!(iline >> ras >> decs))
		{
		  cerr << " ERROR : cannot decode 2 coordinates in :" << line;
		  status = EXIT_FAILURE;
		  break;
		}
	      double ra = RaStringToDeg(ras); // in principe decodes both sexagesimal and decimal.
	      double dec = DecStringToDeg(decs); // same comment.
	      double x,y;
	      RaDec2Pix->apply(ra, dec, x, y);
	      if (skip_out_of_bounds && !frame.InFrame(x,y)) continue;
	      cout << x +MEMPIX2DISK << " " << y + MEMPIX2DISK << endl;
	    }
	  fclose(list);
	}
    }
  return status;
}
