#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "dicstar.h"
#include "listmatch.h"
#include "gtransfo.h"
#include "starmatch.h"

static void usage(const char *prog)
{
  std::cerr << prog << " list1 list2 -o <outname> -c <cut>" << std::endl
	    << " where cut is to be given in units of list coordinates" 
	    << std::endl;
  exit(1);
}



int main(int nargs, char **args)
{
  std::vector<std::string> names;
  double matchCut = 0;
  std::string outName;
  for (int i=1; i <nargs; ++i)
    {
      char *arg = args[i];
      if (arg[0] != '-')
	{
	  names.push_back(arg);
	  continue;
	}
      switch (arg[1]) 
	{
	case 'c' : ++i; matchCut = atof(args[i]); break;
	case 'o' : ++i; outName = args[i]; break;
	default: std::cerr << "don't understand " << arg << std::endl; 
	  usage(args[0]); break;
	}
    }
  if (names.size() > 2) 
    {
      std::cerr << " matching more than 2 lists is not implemented yet" << std::endl;
      exit(1);
    }
  if (names.size() != 2 || outName == "" || matchCut == 0) usage(args[0]);

  GtransfoIdentity id;
  DicStarList l1(names[0]);
  std::cout << " read " << l1.size() << " objects in " << names[0] << std::endl;
  DicStarList l2(names[1]);
  std::cout << " read " << l2.size() << " objects in " << names[1] << std::endl;
  StarMatchList *sm = ListMatchCollect((const BaseStarList &) l1,
				       (const BaseStarList &) l2,
				       &id, matchCut);

  std::cout << " collected " << sm->size() << " matches " << std::endl;
  sm->write(outName);
  delete sm;
  return 1;
}
