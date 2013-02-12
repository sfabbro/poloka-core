#include <iostream>
#include <sstream>
#include <cstdlib> // for strtod
#include <cmath> // for fabs 
#include <iomanip> // for setprecision

#include <poloka/fitsimage.h>
#include <poloka/wcsutils.h>
#include <poloka/astroutils.h>
#include <poloka/fileutils.h>
#include <poloka/gtransfo.h>
#include <poloka/basestar.h>
#include <poloka/frame.h>

static void usage(const char* progname)
{
  cerr << "Usage: " << progname << " [OPTION]... [-c|-p ARGS] FITS...\n"
       << "Convert x-y to sky coordinates and vice-versa using WCS information\n\n"
       << "  -c A D  : convert alpha delta to x-y\n"
       << "  -c @FILE: convert alpha delta to x y from FILE (one alpha delta per line)\n"
       << "  -p X Y  : convert x-y to alpha delta\n"
       << "  -p @FILE: convert x-y to alpha delta from FILE (one x y per line)\n"
       << "  -d      : use decimal degrees\n"
       << "  -s      : skips out of frame coordinates (use with @FILENAME)\n"
       << "  -v      : dump WCS information\n\n";
  exit(EXIT_FAILURE);
}


static void convert_direct(const GtransfoRef Pix2RaDec, 
			   const double& x, const double& y,
			   const bool decimal,
			   const char* Remainder = 0)
{
  double ra,dec;
  Pix2RaDec->apply(x,y,ra,dec);
  if (decimal)
    cout << setprecision(10) << ra << " " << dec;
  else
    cout << RaDegToString(ra) << " " << DecDegToString(dec);
  if (Remainder)
    cout << " " << Remainder;
  else 
    cout << endl;
}

struct Coordinates {
  double x,y;
  string ras, decs;

  // put a default constructor to enventually trace uninitialized data
  Coordinates() : x(-1), y(-1), ras("-1"), decs("92") {}
};


int main(int argc, char **argv)
{
  if (argc < 4) usage(argv[0]);

  list<string> fitsList;
  double x;
  double y;
  bool pix2angles = false;
  bool angles2pix = false;
  bool decimal = false;
  bool skip_out_of_bounds = false;
  bool verbose = false;
  vector<Coordinates> to_convert;
  FILE *listFile = NULL;
  int status = EXIT_SUCCESS;

  for (int i=1; i< argc; ++i)
    {
      char *arg = argv[i];
      
      if (arg[0] != '-')
	{
	  fitsList.push_back(arg);
	  continue;
	}
      switch (arg[1])
	{
	case 'p' :
	  pix2angles = true;
	  i++; 
	  if (argv[i][0] =='@')
	    {
	      const char* fileName = argv[i]+1; 
	      listFile = fopen(fileName,"r");
	      if (!listFile) 
		{
		  cerr << argv[0] << ": cannot open " << fileName << endl; 
		  return EXIT_FAILURE;
		}
	      break;
	    }
	  else 
	    {
	      if (angles2pix)
		{
		  cerr << argv[0] << ": cannot mix convertions both ways, choose -p or -c\n";
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
		  listFile = fopen(fileName,"r");
		  if (!listFile) 
		    {
		      cerr << argv[0] << ": cannot open " << fileName << endl;
		      return EXIT_FAILURE;
		    }
		  break;
		}
	      else
		if (pix2angles)
		  {
		    cerr << argv[0] << ": cannot mix convertions both ways, choose -p or -c\n";
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
	  cerr << argv[0] << ": bad argument " << argv[i] << endl;
	  usage(argv[0]);
	  return EXIT_FAILURE;
	}
    } // end looping on args

  for (list<string>::const_iterator it = fitsList.begin(); it != fitsList.end(); ++it) {
    string fitsName = *it;
    if (!FileExists(fitsName))
      {
	cerr << argv[0] << ": cannot open " << fitsName << endl;
	return EXIT_FAILURE;
      }
    FitsHeader head(fitsName);
    GtransfoRef Pix2RaDec = WCSFromHeader(head);
    if (!Pix2RaDec)
      {
	cerr << argv[0] << ": do not find the expected WCS in " << fitsName  << endl;
	return 0;
      }
    if (verbose) cout << " WCS : " << endl << *Pix2RaDec << " ======= " << endl; 
    
    
    if (pix2angles)
      {
	if (!listFile) 
	  {
	    for (size_t k=0; k < to_convert.size(); ++k)
	      {
		Coordinates &c = to_convert[k]; 
		convert_direct(Pix2RaDec, c.x, c.y, decimal);
	    }
	  }
	else
	  {
	    char line[8192];
	    char* remainder = NULL;
	    while(fgets(line, 8192, listFile))
	      {
		if (line[0] == '#')
		  {
		    if (strstr(line,"format"))
		      {
			int format = DecodeFormat(line,"BaseStar");
			if (format != 2) {
			  cerr << argv[0] << ": format " << format << " invalid\n";
			  return EXIT_FAILURE;
			}
		      }
		    cout << line; continue;
		  }
		if (line[0] == '@') {cout << line; continue;}
		char *next_stuff = line;
		x = strtod(line, &next_stuff);
		if (next_stuff != line) y = strtod(next_stuff, &remainder);
		if (next_stuff == remainder)
		  {
		    cerr << argv[0] << ": cannot decode 2 coordinates in :" << line;
		    status = EXIT_FAILURE;
		    break;
		  }
		convert_direct(Pix2RaDec,x,y,decimal, remainder);
	      }
	    fclose(listFile);
	  }
      }
    else  if (angles2pix)
      {
	// 0.01 is the precision ot invertion (in pixels);
	Frame frame(head);
	GtransfoRef RaDec2Pix = Pix2RaDec->InverseTransfo(0.01, frame);
	if (!listFile) 
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
		    cerr << argv[0] << ": too large error when inverting WCS\n";
		    status = EXIT_FAILURE;
		  }
	      }
	  }
	else
	  {
	    char line[8192];
	    string ras,decs;
	    while(fgets(line, 8192, listFile))
	      {
		if (line[0] == '#' || line[0] == '@') continue;
		istringstream iline(line);	      
		if (!(iline >> ras >> decs))
		  {
		    cerr << argv[0] << ": cannot decode 2 coordinates in :" << line;
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
	    fclose(listFile);
	  }
      }
  }
  return status;
}
