#include "reducedimage.h"
#include "daoutils.h"

static void usage(char *progName)
{
  cerr << progName << " <DbImage>  :  builds files to run DAOPHOT/ALLSTAR" << endl;
}


int main(int nargs, char **args)
{
  if (nargs < 2)  {usage(args[0]);  exit(1);}
  bool overwrite = false;
  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      if ((arg[0] == '-') && (arg[1] == 'o')) overwrite = true;  
      else 
	{
	  string name = string(args[i]); 
	  ReducedImage dbim(name);
	  if (!Sex2Dao(dbim,overwrite)) {cerr << " sex2dao " << name << " : FAILURE!! "<< endl;}
	}
    }
  return EXIT_SUCCESS;
}
