#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

#include <poloka/dicstar.h>
#include <poloka/listmatch.h>
#include <poloka/gtransfo.h>
#include <poloka/nstarmatch.h>
#include <poloka/polokaexception.h>

static void usage(const char *progname)
{
  cerr << "Usage :" << progname << "[OPTION]... LIST...\n"
       << "Join catalogues into one, optionnally filtering entries\n\n"
       << "    -c FLOAT : cut is to be given in units of list coordinate's units\n"
       << "    -t STRING: tags characters are used to tag output in the order of input lists\n"
       << "    -n INT   : min number of matches (default is number of lists)\n"
       << "    -a INT   : ambiguities_handling can be 0,1,2 or 3 (default is 3)\n"
       << "    -f       : will also find transformation between lists\n"
       << "    -o FILE  : output list to FILE (default is joined.list)\n\n";
  exit(EXIT_FAILURE);
}


int main(int nargs, char **args) {

  if (nargs < 3) usage(args[0]);

  vector<string> imList;
  double matchCut = 0;
  string outName = "joined.list";
  string tags;
  int minCount = 0;
  bool fittransfo = false;
  int remove_ambiguities = 3;

  // args collecting
  for (int i=1; i <nargs; ++i) {
    char *arg = args[i];
    if (arg[0] != '-') {
      imList.push_back(arg);
      continue;
    }
    switch (arg[1]) {
    case 'c' : ++i; matchCut = atof(args[i]); break;
    case 'a' : ++i; remove_ambiguities = atoi(args[i]);break;
    case 'o' : ++i; outName = args[i]; break;
    case 't' : ++i; tags = args[i]; break;
    case 'n' : ++i; minCount = atoi(args[i]); break;
    case 'f' : fittransfo = true; break;
    default: cerr << "don't understand " << arg << endl; 
      usage(args[0]); break;
    }
  }

  // args checking
  if (outName.empty()) {
    cerr << args[0] << ": need an output file name\n";
    return EXIT_FAILURE;
  }

  if (matchCut == 0) {
    cerr << args[0] << ": need a matchcut value\n";
    return EXIT_FAILURE;
  }

  if (imList.size() <= 1) {
    cerr << args[0] << ": need t least 2 input lists\n";
    return EXIT_FAILURE;
  }

  if (!tags.empty() && tags.size() < imList.size()) {
    cerr << args[0] << tags.size() << " tags should be >= " << imList.size() << endl;
    return EXIT_FAILURE;
  }


  /* if 2 input lists are provided restore the "old" 2 lists
     behaviour (when joinlist only handled 2 input lists) */
  if (minCount == 0 /* not provided */ && imList.size() == 2)
    minCount = 2;
  
  try {
    NStarMatchList nsm;
    for (size_t k=0; k<imList.size(); ++k) {
      DicStarList l(imList[k]);
      if (l.empty()) {
	cerr << args[0] << ": skipping " << imList[k] << " because it looks empty\n";
	continue;
      }
      cout << " read " << l.size() << " objects from " << imList[k] << endl;
      string tag;
      if (!tags.empty()) tag = tags.substr(k,1);
      if (fittransfo && k > 0) {
	GtransfoRef transfo = ListMatch((BaseStarList&)l, (BaseStarList&)nsm);
	cout << " found transfo " << *transfo << endl;
	TStarList tlist((const BaseStarList &) l, *transfo);
	TStar *empty = new TStar(*(l.EmptyStar()), *transfo);
	nsm.MatchAnotherList((const BaseStarList &) tlist, matchCut, 
			     empty, tag);	
      } else {
	nsm.MatchAnotherList((const BaseStarList &) l, matchCut,
			     l.EmptyStar(), tag, remove_ambiguities);
      }
    }
  
    if (minCount) nsm.ApplyCountCut(minCount);
  
    cout << " writing " << nsm.size() << " matches to " << outName << endl;
    nsm.write(outName);
  } catch(PolokaException p) {
    p.PrintMessage(cout);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
