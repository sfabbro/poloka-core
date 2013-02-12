#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include <poloka/seeing_box.h>
#include <poloka/sextractor_box.h>
#include <poloka/reducedimage.h>
#include <poloka/fileutils.h>
#include <poloka/polokaexception.h>



static void usage(const char *progName)
{
  cerr << "Usage: " << progName << " [OPTION]... [DBIMAGE]...\n"
       << "Perform SExtractor cataloguing and seeing computation\n\n"
       << "  -O : redo from beginning, re-adding background\n"
       << "  -o : overwrite\n"
       << "  -S : write saturated pixel map\n"
       << "  -N : do not subtract background whether it is used or not\n"
       << "  -d : use the sigma background in header for detection\n\n";
  exit(EXIT_FAILURE);
}


int main(int argc, char **argv)
{
 
  if (argc<2) usage(argv[0]);

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
	      usage(argv[0]);
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
	      usage(argv[0]);
	    }
	}
    }/* end args loop */



 
  
  DbImageList list(names); // expands correctly...
  int status = 1;
  
  for (DbImageIterator it = list.begin(); it != list.end(); ++it)
    {

      try { 
	
	ReducedImage redimage(*it);
	if (specif)
	  {
	    status &= redimage.MakeCatalog(redo_from_beg, overwrite, savemasksat,pas_sub_fond,
					   use_sigma_header);
	    redimage.MakeCosmic();
	  }
	else
	  status &= redimage.MakeCatalog();

      }catch(PolokaException p) {
	p.PrintMessage(cout);
	status = 0;
      }
    }
  if (status) return EXIT_SUCCESS; else return EXIT_FAILURE;

}
