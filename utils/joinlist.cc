#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "dicstar.h"
#include "listmatch.h"
#include "gtransfo.h"
#include "nstarmatch.h"

static void usage(const char *prog)
{
  cerr << " usage :" << endl;
  cerr << prog << " list1 list2 ... listn -o <outname> -c <cut> [-t <tags>] [-n <minMatches>]" << std::endl
       << " where : " << endl
       <<" - cut is to be given in units of list coordinate's units" 
       << std::endl
       << " - tags characters are used to tag output in the order of input lists (eg -t griz) " << endl
       << " - matches with less than <minCount> are dropped if provided" 
       << "     (set to 2 for 2 input lists by default)" << endl;
  
  exit(1);
}


int main(int nargs, char **args)
{
  std::vector<std::string> names;
  double matchCut = 0;
  std::string outName;
  std::string tags;
  int minCount = 0;
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
	case 't' : ++i; tags = args[i]; break;
	case 'n' : ++i; minCount = atoi(args[i]); break;
	default: std::cerr << "don't understand " << arg << std::endl; 
	  usage(args[0]); break;
	}
    }

  // args checking
  if (outName == "")
    {
      cerr << " need an output file name !! " << endl;
      usage(args[0]);
      }

  if (matchCut == 0)
    {
      cerr << " need a matchcut value " << endl;
      usage(args[0]);
    }

  if (names.size() <= 1)
    {
      cerr << " you have to provide at least 2 input lists !!! " << endl;
      usage(args[0]);
    }

  if (tags != "" && tags.size() < names.size())
    {
      cerr << " ERROR : you should provide at least as many tags as lists " << endl
	   << " you provided " << tags.size() << " tags '" << tags << "' " 
	   << " for " << names.size() << " input lists " << endl;
      usage(args[0]);
    }
  // arguments are hopefully correct      


  /* if 2 input lists are provided restore the "old" 2 lists
     behaviour (when joinlist only handled 2 input lists) */
  if (minCount == 0 /* not provided */ && names.size() == 2)
    minCount = 2;


  NStarMatchList nsm;
  for (unsigned k=0; k < names.size(); ++k)
    {
      DicStarList l(names[k]);
      if (l.empty())
	{
	  cerr << " skipping " << names[k] << " because it looks empty " << endl;
	  continue;
	}
      cout << " read " << l.size() << " objects from " << names[k] << endl;
      string tag;
      if (tags != "") tag = tags.substr(k,1);
      nsm.MatchAnotherList((const BaseStarList &) l, matchCut, 
			   l.EmptyStar(), tag);
      //      if (k==1) check_list(nsm, " dans la boucle");
    }

  if (minCount) nsm.ApplyCountCut(minCount);


  //  check_list(nsm, " sortie de la boucles " );

  cout << " writing " << nsm.size() << " matches to " << outName << endl;
  nsm.write(outName);

#ifdef STORAGE


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
#endif


  return EXIT_SUCCESS;
}
