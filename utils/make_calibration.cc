#include "reducedimage.h"
#include "jimstar.h"
#include "fidstar.h"
#include "photstar.h"
#include "lightcurvebuilder.h"
#include "listmatch.h"
#include "starmatch.h"


//typedef StarList<LightCurve> LightCurveList;

/*! Dump help for this program */
void usage(const char* programname) {
  
  cout << " syntax  : " << programname << " <refstars.list> <reference reduced image> <catalog>" << endl;
  exit(-1);
}


int main(int argc, char **argv)
{
  if (argc < 4) {
    usage(argv[0]);
  }
  
  string refstarsfile  = argv[1];
  string redimagename = argv[2];
  string catalog = argv[3];
  
  
  ReducedImage redimage(redimagename);
  FitsHeader header(redimage.FitsImageName(Calibrated));
  
  JimStarList* jimlist = GetSelectedJimStarList(header,redimage.UsablePart(),catalog);
  // now x and y of jimstars are already in the reference image frame
  
  // let's read the catalogue of reference stars 
  FidStarList reflist(refstarsfile);
  cout << "reflist.size()=" << reflist.size() << endl;
  
  FidStarCIterator si = reflist.begin();
  for(unsigned int i=0;i<reflist.size(); i++) {
    cout << (**si) << endl;
    ++si;
  }
  
  // ok now we match both catalogs
  GtransfoIdentity id;
  StarMatchList *matchlist = ListMatchCollect((BaseStarList&)(*jimlist),(BaseStarList&)reflist,&id,5);
  cout << "Size of matchlist = " << matchlist->size() << endl;
  string matchlistname = refstarsfile;
  matchlistname += ".match";
  matchlist->write(matchlistname);

  return 0;
}

