#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <string>


#include "seeing_box.h"
#include "dbimage.h"
#include "reducedimage.h"
#include "fileutils.h"



void usage(const char* pg)
{
  cerr << "usage: " << pg << " <dbimage1> <dbimage2> ... ";
  cerr << "-o: overwrite fits key SESEEING" << endl ;
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
 
  if(argc<2)
    usage(argv[0]);
  bool overwrite = false ;
  list<string> names;
  for (int i=1; i< argc; i++) {
    char *arg = argv[i];
    if (arg[0] != '-') {
      names.push_back(arg);
    } else {
      switch (arg[1]) {
      case 'o' :
	overwrite = true;
	break;
      default:
	cerr << "Unknown option " << arg << endl;
	usage(argv[0]);
      }
    }
  }
  DbImageList list(names);
  for (DbImageIterator it = list.begin(); it != list.end(); ++it) {
    ReducedImage redimage(*it);
    if(! redimage.HasCatalog() ) {
      cerr << "dbimage " << redimage.Name() << " has no catalog" << endl;
      return EXIT_FAILURE;
    }
    SEStarList stlse(redimage.CatalogName());
    DatSeeing datsee(redimage.Saturation()) ;
    SortieSeeing sortiese;
    SEStarList seestar;
    CalculeSeeingSE(datsee, sortiese, stlse, seestar);
    if(overwrite) {
      cout << redimage.Name() << " write_seeing " << sortiese.seeing << endl;
      redimage.SetSeeing(sortiese.seeing,"Poloka Weighted Histogram Method");
    }else{
      float oldseeing = redimage.Seeing();
      cout << redimage.Name() << " old_new " << oldseeing
	   << " " << sortiese.seeing;
      if(oldseeing>0)
	cout << " " << sortiese.seeing/oldseeing-1 << endl;
      else
	cout << endl;
    }
  }
  return EXIT_SUCCESS;
}

