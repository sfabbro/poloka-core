#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "newsub.h"
#include "fileutils.h"

#ifdef __sun
#define __GNU_LIBRARY__
/* #define __STDC__ */
/*
extern "C" {
extern int getopt(int, char *const *, const char *);
}
*/
/* #define __EXTENSIONS__ deos not work ?? */
#endif
#include <getopt.h>

void usage()
{
  cerr << "usage: 1 sub at a time for the moment !!!!" ; 
  cerr << "-h: usage" << endl ;
  cerr << "-o: overwrite" << endl ;
  cerr << "-d: detection only (all output images are present)" << endl ;
  cerr << "-c: re-compute parameters and do cuts on already detected candidates only (output images needed)" << endl ;  
  cerr << "-C: cuts on already detected candidates only (no images needed)" << endl ;
  cerr << "-S: images wedding file" << endl ;
  cerr << "-D: sub datacard (default - TOADSCARDS/sub.datacard)" << endl ;
  cerr << "-P: run Pierre's subtraction " << endl;
  cerr << "-x  -y: coordinate of candidate to look for " << endl;
  cerr << " ca m'etonnerais que tout ca marche ! (P.A) " << endl;
  exit(0);
}


static void 
SubtractionProcess(const string fichier, int overwrite, const bool OnlyDet)
{


  NewSub sub(fichier, overwrite, OnlyDet);  
  //  sub.Set_Candidate_Coordinates(xcand,ycand ) ;

#ifdef STORAGE
  if ( ( overwrite == 0 ) && ( FileExists(sub.CandCut_List()) ) )
    {
      cerr << " candidate catalog " << sub.CandCut_List() << " already done for subtraction file  " << fichier 
	   << endl ;
      return;
    }

  if ( only_det  )
    { 
      cerr << " doing only the candidates detection " << endl ;
      sub.OnlyDetection();
      return;
    }

  if ( only_constr  )
    { 
      cerr << " computing parameters and doing cuts on already found candidates  " << endl ;
      sub.OnlyConstruction();
      return;
    }

  if ( only_cut  )
    { 
      cerr << " doing only thecuts on already found candidates  " << endl ;
      sub.OnlyCut();
      return;
    }


  if (p_ou_d == Pierre)   sub.DoItPierre();
  if (p_ou_d == Delphine ) sub.DoIt();
#endif
  sub.DoIt();
  return;

}

int
main(int argc, char **argv)
{
  char c;

  string fichier = "subfile";
  string nomdat;
  int overwrite = 0 ;
  int ouidat = 0 ;
  bool only_det = false , only_constr = false, only_cut = false;

  double xcand = -1, ycand = -1 ;

  while ((c = getopt(argc, argv, "hS:oD:dx:y:Cc")) != -1) 
    {
      switch (c)
	{
	case 'h' :
	  usage();
	  break;
	case 'o' :
	  overwrite =1;
	  break;
	case 'x' :
	  sscanf(optarg, "%lf", &xcand);
	  break;
	case 'y' :
	  sscanf(optarg, "%lf", &ycand);
	  break;
	case 'd' :
	  only_det = true;
	  break;
	case 'c' :
	  only_constr = true;
	  break;
	case 'C' :
	  only_cut = true;
	  break;
	case 'S' :
	  fichier = optarg ;
	  break;
	case 'D' :
	  nomdat = optarg ; ouidat = 1 ;
	  break;
	default:
	  usage();
	}
    }
    
  if (optind >= argc+1) usage();

  if ( ouidat == 0) // non specified datacard
    {
      cerr << "datacard name non specified. " << endl ;
      char *datacard_dir = getenv("TOADSCARDS");
      if (datacard_dir == NULL )
	{
	  cerr << "TOADSCARDS not defined, wont use any datacard  " 
	       << endl;
	}
      else
	{
	  string sdir = datacard_dir ;
	  string nomdatacard = sdir + "/sub.datacard" ;
	  if ( !FileExists(nomdatacard.c_str()) )
	    cerr << "datacard:" 
		 << nomdatacard << " not found, wont use any datacard  " 
		 << endl;
	  else
	    {
	      nomdat =  nomdatacard ;
	      cerr << "will use : " << nomdat << endl ;
	    }

	}
    }
  SubtractionProcess(fichier, overwrite, only_det);
}
