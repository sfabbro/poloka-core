#include <sstream>
#include <iostream>
#include <vector>

#include "stringlist.h"
#include "frame.h"
#include "gtransfo.h"
#include "fileutils.h"
#include "astroutils.h"
#include "wcsutils.h"
#include "fitsimage.h"
#include "dbimage.h"
#include "imageutils.h"

static void usage(const char *progName) {
  cerr << "Usage: " << progName << " [OPTION]... FITS...\n"
       << "Select and list FITS images with criteria of keys in their header\n\n"
       << "   -cal           : use dbimages (calibrated, raw)\n"
       << "   -date [MIN,MAX]: UTC date interval, format is YYYY-MM-DD or MJD\n"
       << "   -band C        : select by band (given by TOADBAND)\n"
       << "   -chip N        : select by chip number (given by TOADCHIP)\n"
       << "   -ra RA         : select by RA (HH:MM:SS or degrees)\n"
       << "   -dec DEC       : select by DEC (DD:MM:SS or degrees)\n"
       << "   -rad VAL       : specify radius (arcmin) from RA and DEC (default: 1')\n"
       << "   -area VAL      : min area (arcmin2) to match image set by -ref (default: 20)\n"
       << "   -ref FITS      : select such that WCS matches a reference FITS\n"
       << "   -s 'KEY OP VAL': select from a general expression:\n"
       << "            KEY -> FITS key to use\n"
       << "            OP  -> operator (==,<,>,<=,>=,~=) ~= is search string\n"
       << "            VAL -> FITS value to match\n";
  exit(EXIT_FAILURE);
}

bool DecodeFitsExpression(const FitsHeader &Head, const string&Exp) {
  istringstream iss(Exp);
  string keyname, op, keyval;
  iss >> keyname >> op >> keyval;
  if (!Head.HasKey(keyname)) return false;
  if (op == "==") 
    return string(Head.KeyVal(keyname)) == keyval;
  if (op == "~=")
    return string(Head.KeyVal(keyname)).find(keyval) != string::npos;
  double dval = atof(keyval.c_str());
  if (op == "<") 
    return double(Head.KeyVal(keyname)) <  dval;
  if (op == "<=") 
    return double(Head.KeyVal(keyname)) <= dval;
  if (op == ">") 
    return double(Head.KeyVal(keyname)) >  dval;
  if (op == ">=") 
    return double(Head.KeyVal(keyname)) >= dval;
  return true;
}

struct ImageSelector {

  double mjdMin, mjdMax, minArea;
  StringList expList;
  bool withDbImages;  
  Frame *radecFrame;
  FitsHeader *refHead;

  ImageSelector() 
    : mjdMin(-1), minArea(20), withDbImages(false),
      refHead(0), radecFrame(0) {}

  ~ImageSelector() {
    if (refHead) delete refHead;
    if (radecFrame) delete radecFrame;
  }

  void operator () (const string& imname) const {
    string fitsname = imname;
    // get the FITS file
    if (withDbImages) {
      if (!DbImage(imname).IsValid()) return;      
      fitsname = DbImage(imname).FitsImageName(Calibrated);
      if (!FileExists(fitsname)) fitsname = DbImage(imname).FitsImageName(Raw);
    }
    if (!FileExists(fitsname)) return;

    FitsHeader head(fitsname);
    if (!head.IsValid()) return;

    // decode dates
    bool filter_date = true;
    if (mjdMin>0) {
      double mjd = ModJulianDay(head);
      filter_date = (mjd >= mjdMin && mjd <= mjdMax);
    }
      
    // decode ref matching
    bool filter_ref = true;
    if (refHead)
      filter_ref = (Arcmin2Overlap(*refHead, head) >= minArea);

    // decode radius
    bool filter_rad = true;
    if (radecFrame) {
      Frame pixFrame(head);
      GtransfoRef Pix2RaDec = WCSFromHeader(head);
      filter_rad = false;
      if (Pix2RaDec) {
	Frame overlap =  ApplyTransfo(pixFrame, *Pix2RaDec) * (*radecFrame);
	filter_rad = (overlap.Area()>0);
      } else
	cerr << head.FileName() << ": missing WCS\n";
    }
    
    // decode expressions
    bool filter_exp = true;
    for (StringCIterator it=expList.begin(); it != expList.end(); ++it)
      filter_exp &= DecodeFitsExpression(head, *it);
    
    // print image if selected
    if (filter_exp && filter_date && filter_rad && filter_ref)
      cout << imname << endl;
  }

};

int main(int nargs, char **args) {

  if (nargs<2) usage(args[0]);

  string datestr(""), refstr("");
  string rastr(""), decstr("");
  StringList toSelect;
  double radius = 1;
  ImageSelector imSelect;

  // loop over arguments
  for (int i=1; i<nargs; ++i) {
    char *arg = args[i];

    // options
    // images
    if (arg[0] == '-') {
      arg++;
      if (strcmp(arg,"cal") == 0) { 
	imSelect.withDbImages = true; 
	continue; 
      }
      if (strcmp(arg,"band") == 0) {
	imSelect.expList.push_back(string("TOADBAND == ")+args[++i]);
	continue; 
      }
      if (strcmp(arg,"chip") == 0) {
	imSelect.expList.push_back(string("TOADCHIP == ")+args[++i]);
	continue;
      }
      if (strcmp(arg,"s") == 0) {
	imSelect.expList.push_back(args[++i]);
	continue;
      }
      if (strcmp(arg,"date") == 0) {
	datestr = args[++i];
	continue;
      }
      if (strcmp(arg,"ref") == 0) {
	refstr = args[++i];
	continue;
      }
      if (strcmp(arg,"area") == 0) {
	imSelect.minArea = atof(args[++i]);
	continue;
      }
      if (strcmp(arg,"rad") == 0) {
	radius = atof(args[++i])/60./* convert to deg */;
	continue; 
      }
      if (strcmp(arg,"ra") == 0) {
	rastr  = args[++i];
	continue;
      }
      if (strcmp(arg,"dec") == 0) {
	decstr = args[++i];
	continue;
      }
    } else {
      toSelect.push_back(arg);
      continue;
    }
    
    // unrecognized option
    usage(args[0]);      
  }

  // decode dates
  if (!datestr.empty()) {	
    int d1,m1,y1,d2,m2,y2;
    if (sscanf(datestr.c_str(),"[%d-%d-%d,%d-%d-%d]", &y1,&m1,&d1,&y2,&m2,&d2) == 6) {
      imSelect.mjdMin = JulianDay(d1,m1,y1) - 2400000.5;
      imSelect.mjdMax = JulianDay(d2,m2,y2) - 2400000.5;
    } else if (sscanf(datestr.c_str(),"[%f,%f]", &imSelect.mjdMin, &imSelect.mjdMax) != 2) {
      cerr << datestr << ": wrong date format\n";
      return(EXIT_FAILURE);
    }
  }
  
  // decode reference
  if (!refstr.empty()) {
    if (imSelect.withDbImages) {
      refstr = DbImage(refstr).FitsImageName(Calibrated);
      if (!FileExists(refstr)) refstr = DbImage(refstr).FitsImageName(Raw);
      if (!FileExists(refstr)) {
	cerr << refstr << ": not a valid reference file\n";
	return(EXIT_FAILURE);
      }
    }
    imSelect.refHead = new FitsHeader(refstr);
  }

  // decode ra, dec, radius
  if (!rastr.empty() && !decstr.empty()) {
    double ra;
    if (sscanf(rastr.c_str(),"%f", &ra) != 1) {
      ra = RaStringToDeg(rastr);
      if (ra >= 99) {
	cerr << rastr << ": wrong RA format\n";
	return(EXIT_FAILURE);
      }
    }
    double dec;
    if (sscanf(decstr.c_str(),"%f", &dec) != 1) {
      dec = DecStringToDeg(rastr);
      if (dec >= 199) {
	cerr << decstr << ": wrong DEC format\n";
	return(EXIT_FAILURE);
      }
    }
    double cosdec = cos(M_PI*dec/180);
    imSelect.radecFrame = new Frame(ra - radius*cosdec,
				    dec - radius,
				    ra + radius*cosdec,
				    dec + radius);
  }

  for_each(toSelect.begin(), toSelect.end(), imSelect);

  return EXIT_SUCCESS;
}
