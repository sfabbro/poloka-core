#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <string>


#include "seeing_box.h"
#include "sextractor_box.h"
#include "dbimage.h"
#include "reducedimage.h"
#include "fileutils.h"



void usage()
{
  cerr << "usage: " ; 
  cerr << "-O: redo from beginning, re-adding background." << endl ;
  cerr << "-o: overwrite" << endl ;
  cerr << "-S: write saturated pixel map" << endl ;
  cerr << "-N: background will NOT be subtracted, wether it is used or not " << endl ;
  cerr << "-d use the sigma background in header for detection " << endl ;
  cerr << " <DbImages ... > " << endl;
  cerr << "What is done is:" << endl ;
  cerr << " sextractor catalog and seeing computation if not already done,except when -o used" << endl ;
  cerr << " Datacards are read in $TOADSCARDS " << endl;
  exit(-1);
}




int main(int argc, char **argv)
{
 
  bool overwrite = false ;
  bool savemasksat = false ;
  bool use_sigma_header = false ;
  bool pas_sub_fond = false ;
  bool specif = false ;
  bool redo_from_beg  = false ;
  list<string> names;

  for (int i=1; i< argc; i++)
    {
      char *arg = argv[i];
      if (arg[0] != '-') 
	{
	  names.push_back(arg);
	}
      else
	{
	  switch (arg[1])
	    {
	    case 'h' :
	      usage();
	      break;
	    case 'd' :
	      use_sigma_header = true ;specif = true ;
	      break;
	    case 'o' :
	      overwrite = true ;specif = true ;
	      break;
	    case 'O' :
	      redo_from_beg= true ;specif = true ;
	      break;
	    case 'S' :
	      savemasksat = true ;specif = true ;
	      break;
	    case 'N' :
	      pas_sub_fond = true ; specif = true ;
	      break;
	    case 't' :
	      cerr << " argument -t to make_catalog obsoleted " << endl;
	      break;
	    default:
	      usage();
	    }
	}
    }/* end args loop */



 
  
  DbImageList list(names); // expands correctly...
  int status = 1;

  for (DbImageIterator it = list.begin(); it != list.end(); ++it)
    {
      ReducedImage redimage(*it);
      if (specif)
	{
	  status &= redimage.MakeCatalog(redo_from_beg, overwrite, savemasksat,pas_sub_fond,
				     use_sigma_header);
	  redimage.MakeCosmic();
	}
      else
	status &= redimage.MakeCatalog();
    }
  if (status) return EXIT_SUCCESS; else return EXIT_FAILURE;

}
