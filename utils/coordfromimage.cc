#include <iostream>
#include <reducedimage.h>
#include <imagematch.h>
#include <gtransfo.h>
#include <dicstar.h>
#include <fitsimage.h>


static void usage(const char* pg) {
  cerr << pg << " coordx coordy dbimage_from dbimage_to" << endl;
  cerr << "option:  -c dicstarlist_from dicstarlist_to : uses these catalogs instead of sextractor ones" << endl;
  exit(1);
}

// we have to keep only entries with different coordinates
static BaseStarList* CalibrationList2BaseStarList(const string& filename) {
  DicStarList dicstarlist(filename);
  BaseStarList* baselist = new BaseStarList();
  Point sp(-1000,-1000);
  for(DicStarCIterator it=dicstarlist.begin(); it !=dicstarlist.end(); ++it) {
    if(sp.Distance(**it)>0.001) {
      CountedRef<BaseStar> star((const BaseStar*)(*it));
      sp = *star;
      baselist->push_back(star);
    }
  }
  return baselist;
}

int main (int argc, char **argv) {
  if (argc < 5)  usage(argv[0]);
  int i=1;
  double coord_x = atof(argv[i++]);
  double coord_y = atof(argv[i++]);
  ReducedImage image_from(argv[i++]);
  ReducedImage image_to(argv[i++]);

  CountedRef<Gtransfo> direct,reverse;
  
  bool specifiedcatalogs=false;
  string catalog_from = "";
  string catalog_to   = "";
  
  if(argc>5) {
  const char* arg = argv[i];
  switch(arg[1]) {
  case 'c':
    catalog_from = argv[++i];
    catalog_to   = argv[++i];
    specifiedcatalogs = true;
    break;
  default:
    cerr << "unexpected argument " << arg << endl;
    usage(argv[0]);
  }
  }
  if(specifiedcatalogs) {
    BaseStarList* list_from =  CalibrationList2BaseStarList(catalog_from);
    BaseStarList* list_to =  CalibrationList2BaseStarList(catalog_to);
    
    FitsHeader head_from(image_from.FitsName());
    FitsHeader head_to(image_to.FitsName());
    
    MatchGuess(*list_from,*list_to,head_from,head_to,direct,reverse);
    RefineGuess(*list_from,*list_to,direct,reverse);

    delete list_from;
    delete list_to;
  }else{
    ImageListMatch(image_from,image_to,direct, reverse);
  }
  BaseStarList* SE2Base(SEStarList * This);

  Point P_from(coord_x,coord_y);
  Point P_to = direct->apply(P_from);
  
  cout << P_to.x << " " << P_to.y << endl;

  return 0;
}
