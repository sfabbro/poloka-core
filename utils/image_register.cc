
#include <string>
#include <cstdlib> // for the setenv stuff
#include <fstream>
#include <cstring>

#include "imagematch.h"
#include "gtransfo.h"
#include "reducedutils.h"

static void usage(char *progName)
{
  cerr << "Usage: " << progName << " [OPTIONS] <DbImage(s)> " << endl;
  cerr << "   register images relatively to a given geometric reference image \n" 
       << "   OPTIONS:\n"
       << "     -geo <name> : indicate the reference image <name>. Default is first one\n"
       << "     -p : print a quick photometric ratio\n"
       << "     -d : dump the following:\n"
       << "            - the matched star list in <name>.<DbImage>.match.list\n "
       << "            - the transfos (direct and reverse) in <name>.<DbImage>.transfo\n";
  exit(-1);
}

int main(int nargs, char **args)
{

  //default options
  // if nothing is given
  if (nargs < 2)  { usage(args[0]); }
  if (nargs == 2) { cerr << " error: register at least 2 images"; usage(args[0]); }
  bool dump = false, phoratio = false;
  string geoName("NOGEO");
  ReducedImageList toRegister;

  //loop over arguments
  for (int i=1; i<nargs; ++i)
    {
      char *arg = args[i];
      // images
      if (arg[0] != '-') 
	{
	  ReducedImage *im = new ReducedImage(arg);
	  toRegister.push_back(im);
	  continue;
	}
      // options
      arg++;
      if (strcmp(arg,"geo")==0) {++i; geoName = args[i]; continue;}
      if (strcmp(arg,"d")==0) { dump = true; continue;}
      if (strcmp(arg,"p")==0) { phoratio = true; continue;}

      // unrecognized option
      usage(args[0]);
    }

  if (geoName=="NOGEO") geoName = toRegister.front()->Name();

  ReducedImage geoRef(geoName);  
  char *old_dump = getenv("DUMP_MATCH");
  for (ReducedImageCIterator it = toRegister.begin(); it != toRegister.end(); ++it)
    {
      const ReducedImage *current = *it;
      CountedRef<Gtransfo> direct, reverse;
      if (dump) setenv("DUMP_MATCH","YES",1);
      //naming of the list dump is taken care in imagematch.cc
      ImageListMatch(geoRef, *current, direct, reverse);
      if (phoratio)
	{
	  double error;
	  cout << " Photometric ratio: "
	       << QuickPhotomRatio(geoRef, *current, error, direct)
	       << " +/- " << error << endl;
	}
      if (dump)
	{
	  ofstream gout((geoName+"."+current->Name()+".transfo").c_str());
	  gout << "# Transfo " << geoName << " to " << current->Name() << endl;
	  direct->dump(gout);
	  gout << "# Transfo " << current->Name() << " to " << geoName << endl;
	  reverse->dump(gout);
	}   
    }

  if (dump && !old_dump) unsetenv("DUMP_MATCH");

  return 1;
}
