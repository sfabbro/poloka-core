#include <iostream>
#include <fstream>
#include <fileutils.h>
#include <dicstar.h>
#include <basestar.h>
#include <vutils.h>
#include <fastfinder.h>
#include <map>
#include <iomanip>

using namespace std;

static void usage(const char *pgname)
{
  cerr << pgname << " -i <inputdicstarlist> -o <outputdicstarlist> -c <secondarycatalog>" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  
  string inputdicstarlistname = "";
  string outputdicstarlistname = "";
  string secondarycatalogname = "";
  
  if (argc < 7)  {usage(argv[0]);}
  for (int i=1; i<argc; ++i)
    {
      char *arg = argv[i];
      if (arg[0] != '-')
	{
	  cerr << "unexpected parameter " << arg << endl;
	  usage(argv[0]);
	}
      switch (arg[1])
	{
	case 'i' : inputdicstarlistname  = argv[++i]; break;
	case 'o' : outputdicstarlistname = argv[++i]; break;
	case 'c' : secondarycatalogname  = argv[++i]; break;
	default : 
	  cerr << "unknown option " << arg << endl;
	  usage(argv[0]);
	}
    }
  
  
  if(!FileExists(inputdicstarlistname)) {
    cerr << "cant find catalog " << inputdicstarlistname << endl;
    usage(argv[0]);
  }
  if(!FileExists(secondarycatalogname)) {
    cerr << "cant find catalog " << secondarycatalogname << endl;
    usage(argv[0]);
  }
  
  DicStarList inputdicstarlist(inputdicstarlistname);
  DicStarList secondarycatalog(secondarycatalogname);
  
  FastFinder finder((const BaseStarList &)secondarycatalog);

  Point point;
  
  double MaxDist = 1./3600.; // 1 arc sec
  int nok = 0;
  int ntot = 0;
  for(DicStarIterator entry=inputdicstarlist.begin();entry!=inputdicstarlist.end();++entry) {
    ntot++;
    point.x =  (*entry)->getval("ra");
    point.y =  (*entry)->getval("dec");
    const BaseStar * bs = finder.FindClosest(point,MaxDist);
    if(!bs) {
      entry = inputdicstarlist.erase(entry);
    }
    nok++;
    const DicStar* calibratedstar = dynamic_cast<const DicStar *>(bs);
    // ok now change the values of the entry
    (*entry)->setval("u",calibratedstar->getval("mu"));
    (*entry)->setval("ue",calibratedstar->getval("emu"));
    (*entry)->setval("g",calibratedstar->getval("mg"));
    (*entry)->setval("ge",calibratedstar->getval("emg"));
    (*entry)->setval("r",calibratedstar->getval("mr"));
    (*entry)->setval("re",calibratedstar->getval("emr"));
    (*entry)->setval("i",calibratedstar->getval("mi"));
    (*entry)->setval("ie",calibratedstar->getval("emi"));
    (*entry)->setval("z",calibratedstar->getval("mz"));
    (*entry)->setval("ze",calibratedstar->getval("emz"));
  }
  
  cout << nok << " matches of " << ntot << " initial stars" << endl;
  inputdicstarlist.write(outputdicstarlistname);
  return EXIT_SUCCESS;
}


