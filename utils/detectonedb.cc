#include <iostream>
#include <stdlib.h>

#include "newsub.h"
#include "toadscards.h"
#include "fitsimage.h"
#include "dbimage.h"
#include "candstar.h"
#include "reducedimage.h"
#include "gtransfo.h"
#include "dodetection.h"


int main(int nargs,char **args)
{ 
 
  string imageName=args[1];
 
  
  
  ReducedImage OneReduced(imageName);
  FitsImage OneFits(OneReduced.FitsName());


  // read the seeing; used during the detection
  double seeing = OneFits.KeyVal("SESEEING");

  
  // use the usable part of the image for the detection
  Frame frame_det = OneReduced.UsablePart();
  Image img = img.Subimage(frame_det);
  


  // build the mask
  //  FitsImage *dead = new FitsImage(OneReduced.FitsDeadName());
  FitsImage *satur= new FitsImage(OneReduced.FitsSaturName());
  
  Image mask = (*satur);
  
  //delete dead;
  delete satur;
  
  
  mask = mask.Subimage(frame_det);
  
  // Read cuts from the datacards
  DatDetec datdet(DefaultDatacards());
  
  
  
  // Detection of the image
  CandidateStarList stl;
  NewCvDetection(img, mask, stl, datdet, seeing, seeing);

  // Detection done within the usable part => shift the catalog
  GtransfoLinShift shift(frame_det.xMin, frame_det.yMin);
  stl.ApplyTransfo(shift);

  // Write the catalog
  string Name = OneReduced.Name() + "/allcand.list";
  stl.write(Name);

  cout << stl.size() << "  candidats selectionnes " << endl ;
}
