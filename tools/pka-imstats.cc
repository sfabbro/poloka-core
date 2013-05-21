#include <iostream>
#include <iomanip>
#include <string>
#include <list>

#include <poloka/fitsimage.h>
#include <poloka/fileutils.h>
#include <poloka/dbimage.h>

static void usage(const char* progname) {
  cerr << "Usage: " << progname << " [OPTION] FITS...\n"
       << "Compute sky, r.m.s., min and max values of FITS image\n\n"
       << "    -what       : for DBIMAGE uses, where what is raw|cal|dead|weight|...\n\n";
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {

  if (argc < 2) usage(argv[0]);
  list<string> imList;
  string which_info = "cal";
  string type = "IMAGE";

  size_t nchar = 10;

  for (int i=1; i<argc; i++) {
    char *arg = argv[i];
    if (arg[0] != '-') 	{
      if (strlen(arg) > nchar) nchar = strlen(arg);
      imList.push_back(arg);
    } else if (strlen(arg)>2) {  /* search for possible : -raw, -dead, -flat, -cal */
      which_info = arg+1;
      type += "(" + which_info + ")";
    } else {
      cerr << argv[0] << ": " << arg << " unknown option\n";
      usage(argv[0]);
    }
  }

  cout << endl << setiosflags(ios::left)
       << setw(nchar) << type 
       << setiosflags(ios::right)
       << setw(9) << "SKY"
       << setw(7) << "SIG"
       << setw(9) << "MIN"
       << setw(9) << "MAX"
       << setw(10) << "DBSTAT"
       << resetiosflags(ios::right)
       << endl << endl;

  bool ok = true;
  for (list<string>::const_iterator it=imList.begin(); it != imList.end(); ++it) {

    string filename = *it;
    DbImage dbimage(filename);
    if (dbimage.IsValid())
      filename = dbimage.GetFileName(which_info.c_str());

    FitsImage image(filename);
  
    if (!image.IsValid()) {
      cerr << argv[0] << ": " << filename << ": invalid file\n";
      ok = false;
      continue;
    }

    Pixel mean,sigma,minv,maxv;
    image.SkyLevel(&mean, &sigma);
    image.MinMaxValue(&minv, &maxv);
    
    cout << setiosflags(ios::left)
	 << setw(nchar) << *it
	 << setiosflags(ios::right) << setiosflags(ios::fixed)
	 << setw(9) << setprecision(2) << mean << ' '
	 << setw(7) << setprecision(2) << sigma << ' '
	 << setw(9) << setprecision(1) << minv << ' '
	 << setw(9) << setprecision(1) << maxv << ' ';

    if (dbimage.IsValid() && FileExists(dbimage.FitsWeightName())) {
      FitsImage weight(dbimage.FitsWeightName());
      if (FileExists(dbimage.FitsSaturName())) {
	FitsImage satur(dbimage.FitsSaturName());
	weight *= 1 - satur;
      }
      cout << setw(6) << ImageAndWeightError(image, weight);
    } else {
      cout << setw(6) << "none";
    }

    cout << resetiosflags(ios::right)
	 << endl;
  }
  
  return ok? EXIT_SUCCESS : EXIT_FAILURE;
}
