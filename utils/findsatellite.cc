#include "cluster.h"
#include "reducedimage.h"
#include "fitsimage.h"
#include "polokaexception.h"

void usage(char * ProgName)
{
  cerr << ProgName << " <dbim> " << endl;
  cerr << "Produce a binary map of the satellites " << endl;
  exit(1);
}

int main(int nargs, char ** argv)
{
  list<string> names;
  if(nargs ==1) usage(argv[0]);
  for (int i=1; i< nargs; i++)
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
	    default:
	      usage(argv[0]);
	    }
	}
    }/* end args loop */


 DbImageList list(names); // expands correctly...
 bool ok = true;
 
  for (DbImageIterator it = list.begin(); it != list.end(); ++it)
    {

      try{

      ReducedImage redimage(*it);
      
      /*
	ClusterList clustList(redimage);
	clustList.Cut();
	Image mask = clustList.Mask();
	FitsHeader header(redimage.FitsName());
	FitsImage(redimage.FitsSatelliteName(), header,mask);
	*/
      redimage.MakeSatellite();
      }catch(PolokaException p) {
	p.PrintMessage(cout);
	ok = false;
      }

    }
  if(ok)
    return EXIT_SUCCESS;
  return EXIT_FAILURE;
}
