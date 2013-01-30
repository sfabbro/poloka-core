#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "dbimage.h"
#include "reducedimage.h"
#include "fileutils.h"



void usage()
{
  cerr << "usage: make_weight <DbImages ... > " << endl;
  exit(-1);
}




int main(int argc, char **argv)
{ 
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
      status &= redimage.MakeWeight();
    }
  if (status) return EXIT_SUCCESS; else return EXIT_FAILURE;
}
