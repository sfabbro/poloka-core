#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <string>

//#include "vignettes.h"
#include "fileutils.h"
#include "candstar.h"
#include "datacards.h"

 void usage();

void usage()
{
  cerr << "usage: pas de choix des noms " << endl ; 
  cerr << "-h: usage" << endl ;
  cerr << "-d: repertoire" << endl ;
  cerr << "***** to re-do cuts on candstarlist  ******" << endl ;
  cerr << "-i: list to cut" << endl ;
  cerr << "-o: cut list (output)" << endl ;
  cerr << "-D: datacards" << endl ;
  exit(0);
}



int
main(int argc, char **argv)
{
  char c;

  string directo ;
  string nomdat, nomi, nomo, mosa="mosaique_recut";
  int ouicut = 0 , ouidat = 0, ouilist=0;
  while ((c = getopt(argc, argv, "hd:D:i:o:")) != -1) 
    {
      switch (c)
	{
	case 'h' :
	  usage();
	  break;
	case 'd' :
	  directo = optarg ;
	  break;
	case 'D' :
	  nomdat = optarg ; ouidat = 1;
	  break;
	case 'i' :
	  nomi = optarg ; ouilist=1;
	  break;
	case 'o' :
	  nomo = optarg ; ouicut = 1 ;
	  break;
	default:
	  usage();
	}
    }
    
  if (optind >= argc+1) usage();
  string listecand;
  if ( ouicut == 1 )
    {
      cout << nomi << " " << nomo << endl ;
      CandStarList stl(nomi);
      cout << "candidats:" <<  stl.size() << endl ;

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
      
      DatDetec datdet ;
      if ( (nomdat.c_str() != NULL) && 
	   (strlen(nomdat.c_str()) > 1 ) &&
	   (FileExists(nomdat.c_str()))   )
	{
	  DataCards data(nomdat.c_str());
	  datdet.LitDataCard(data);
	}
      else
	{
	  cerr << " taking default configuration" << endl ;
	  datdet.Default();
	}
      datdet.Print();
      CandStarList stlcut;
      stl.Cut(stlcut, datdet);
      string nomo1 = nomo + ".list" ;
      string nomo2 = nomo + ".nice" ;
      stlcut.write(nomo1);
      stlcut.write_nice(nomo2);
      cout << "candidats:" <<  stlcut.size() << endl ;
      BuiltMosaique(directo, nomo1,mosa);

    }
  else
    {
      if ( ouilist == 1)
	BuiltMosaique(directo, nomi,mosa);
      else
	BuiltMosaique(directo,mosa);
    }
	

return(0);

}
