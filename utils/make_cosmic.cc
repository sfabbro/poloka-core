#include "reducedimage.h"
#include "fitsimage.h"
#include "fileutils.h"
#include "sestar.h"
#include "seeing_box.h"



static void usage(const string &exec)
{
  cerr << " usage for " << exec << endl;
  cerr << exec << " <ReducedImages...>" << endl 
       << "    -o : overwrite   " << endl;
  exit(-1);
}

int main(int nargs, char ** args)
{
  ReducedImageList RList;
  if (nargs<2) usage(args[0]);
  bool overwrite = false;

  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      if ((arg[0] != '-')) 
	{
	  ReducedImage *current = new ReducedImage(args[i]);
	  if (FileExists(current->FitsName())) 
	    RList.push_back(current);
	  else cerr << "ReducedImage " << args[i] << " not found "<< endl;
	  continue;
	}
      switch (arg[1])
	{
	case 'o' : overwrite = true; break;
	default : usage(args[0]); exit(1);	  
	}
    }
  
  for(ReducedImageIterator it = RList.begin(); it != RList.end(); ++it)
    {
      ReducedImage *current = *it;
      if (overwrite && current->HasCosmic())
	remove((current->FitsCosmicName()).c_str());
      current->MakeCosmic();

      /*// recompute seeing
	string datacards_dir = getenv("TOADSCARDS");
	if (datacards_dir == "" )
	{
	datacards_dir = "." ;
	cout << " TOADSCARDS not defined, trying to find cards in local 
	directory  " 
	<< endl;
	}
	else cout << " using SExtractor datacards from "  
	<< datacards_dir << endl;
	DatSeeing datSeeing(datacards_dir+"/seeing.datacard",
	current->Saturation() );
	SortieSeeing sortie;
	SEStarList seestar;
	SEStarList selist(current->CatalogName());
	CalculeSeeingSE(datSeeing, sortie, selist, seestar);
	current->SetSeeing(sortie.seeing);
      */
    }
  
  return 0;
}
