// -*- C++ -*-
// 
// file: imsum.cc
//       somme de tous les pixels d'une image
// 
#include <libgen.h>
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "fileutils.h"
#include "fitsimage.h"
#include "dictfile.h"



using namespace std;


void usage()
{
  cerr << "usage: imsum OPTIONS <imagelist>" << endl;
  cerr << "OPTIONS:" << endl;
  cerr << " -I <inputlist>" << endl;
  cerr << " -p just print the sum" << endl;
  cerr << " -o <ntuple>" << endl;
  exit(-1);
}



int main(int argc, char** argv)
{
  string inputlistname;
  string outfile = "sumpixels.ntuple";
  vector<string> imnames;
  
  bool print=false;
  
  int i;
  char c;
  while( (c=getopt(argc, argv, "hI:o:p")) != -1)
    switch(c) {
    case 'I':
      inputlistname = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'p':
      print = true;
      break;
    case 'h':
      usage();
    default:
      usage();
    }
  for(i=optind;i<argc;i++)
    imnames.push_back(argv[i]);
  if(inputlistname == "" && imnames.size() == 0)
    usage();
  
  if(imnames.size() == 1 && print) {
      FitsImage  im(imnames[0]);
      double nx = im.Nx(), ny = im.Ny();
      double npix = nx * ny;
      double sum = im.SumPixels();
      cout << sum << endl;
      return 0;
  }
  
  ofstream ofs(outfile.c_str());
  ofs << "# mmjd :" << endl
      << "# seeing :" << endl
      << "# skysig :" << endl
      << "# ccd :" << endl
      << "# path :" << endl
      << "# band :" << endl
      << "# field : " << endl
      << "# nsat : " << endl
      << "# npix : " << endl
      << "# end" << endl;
  
  if(imnames.size() > 0) {
    vector<string>::iterator I;
    for(I=imnames.begin();I!=imnames.end();I++) {
      FitsImage  im(*I);
      double nx = im.Nx(), ny = im.Ny();
      double npix = nx * ny;
      double sum = im.SumPixels();
      cout << " (*) " << *I << " " << im.SumPixels() << "(" << 100 * sum/npix << "%)" << endl;
    }
  }
  
  DictFile ilst(inputlistname);
  DictFileCIterator I;
  for(I=ilst.begin();I!=ilst.end();I++) {
    string pppp = I->Value("path");
    string path = dirname(const_cast<char*>(pppp.c_str()));
    cout << " (*) processing " << path << endl;
    string saturname = path + string("/satur.fits.gz");
    if(!FileExists(saturname)) {
      continue;
    }
    
    FitsImage im(saturname);
    double nx = im.Nx(), ny = im.Ny();
    double npix = nx * ny;
    double nsat = im.SumPixels();
    ofs << I->Value("mmjd") << " "
	<< I->Value("seeing") << " "
	<< I->Value("skysig") << " "
	<< I->Value("ccd") << " "
	<< I->Value("path") << " "
	<< I->Value("band") << " "
	<< I->Value("field") << " "
	<< nsat << " "
	<< npix << " "
	<< endl;
  }
  
  ofs.close();
  
}
