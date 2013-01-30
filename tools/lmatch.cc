#include <iostream>
#include "senearstar.h"
#include "sestar.h"
#include "listmatch.h"
#include "starmatch.h" 
//#include "subtraction.h"
#include "dbimage.h"
#include "imagematch.h"
#include "reducedutils.h"
#include "starlist.h"

int photom_ratio(int nargs, char **argv)
{
  DbImage i1(argv[1]);
  DbImage i2(argv[2]);

  if (!i1.IsValid() && i2.IsValid())
    {
      cerr <<  " peux pas ouvrir " << argv[1] << " || " << argv[2] << endl;
      return -1;
    }
  //  Gtransfo *direct, *reverse;
  // ImageListMatch(i1, i2, direct, reverse,2);
  double error;
  Gtransfo *toto = new GtransfoIdentity;
  QuickPhotomRatio(i1, i2, error, toto);
  return 1;

}


int image_match(int nargs, char **argv)
{
  DbImage i1(argv[1]);
  DbImage i2(argv[2]);

  if (!i1.IsValid() && i2.IsValid())
    {
      cerr <<  " peux pas ouvrir " << argv[1] << " || " << argv[2] << endl;
      return -1;
    }
  Gtransfo *direct, *reverse;
  ImageListMatch(i1, i2, direct, reverse);
  return 1;
}



int main(int nargs, char **argv)
{

  cout << "Reading lists..." << endl;
  
  SEStarList l1;
  
  SEStarList l2;
  
  l1.read(argv[1]);
  l2.read(argv[2]);
  
  if (l2.size() == 0)
    {
      cout << "CA VA PAS DU TOUT " << l2.size() << endl;
      photom_ratio(nargs,argv);
      //      image_match(nargs, argv);
      exit(0);
    }
  
  double n=3;
  double k=3;
#ifdef STORAGE 
  sscanf(argv[3],"%lf",&n);
  sscanf(argv[4],"%lf",&k);
#endif
  cout << "Searching an initial match..." << endl;
  BaseStarList *bl1 = (BaseStarList*) &l1;
  BaseStarList *bl2 = (BaseStarList*) &l2;
  bl1->write("list1.list");
  bl2->write("list2.list");

  MatchConditions conditions;
  conditions.SizeRatio = 1;
  conditions.DeltaSizeRatio = 0.1;
  conditions.NStarsL1 = 70;
  conditions.NStarsL2 = 70;
  
  StarMatchList *matches = MatchSearchRotShiftFlip(*bl1, *bl2, conditions);
  matches->DumpTransfo();
  cout << "Refining 1st order ..." << endl;  
  cout << "Collecting matches less than "<< n  << " pixels away" << endl;  
  StarMatchList *refined =  ListMatchCollect(*bl1, *bl2, 
					     matches->Transfo(),n);
  
  delete matches;

  refined->RefineTransfo(k);
  refined->DumpTransfo();
  
  refined->write(*(refined->Transfo()),"match12.list");

  return 1;
}
