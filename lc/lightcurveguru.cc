#include "lightcurveguru.h"

static void selectImages(ReducedImageList& Images, ReducedImageList& SelectedImages, Pred Cond)
{
  for (ReducedImageIterator itRed=Images.begin(); itRed != Images.end(); )
    {  
      if (Cond) { SelectedImages.push_back(itRed); itRed = Images.erase(itRed); } 
      else ++itRed;
    }
}

LightCurveGuru::LightCurveGuru(const string& LcFileName) : vector<ReducedImageList(false)>
{
  LightCurveInit lcInit;
  if (!lcInit.read(LcFileName)) return; // or exit(1)
  lcInit.Init();
  *this = lcInit.Images;
}

void LightCurveGuru::MonopolizeCPU()
{  
  cout << " LightCurveGuru::MonopolizeCPU() : Starting light curve production \n";
  
  // align all stars with photometric ratio
  AlignStars align(itSet->PhoRef);
  for_each(itSet->begin(), itSet->end(), align);

  // aperture on subtractions to check variable stars
  ApertureSub aper(itSet->PhoRef);
  for_each(ImageSet.begin(), ImageSet.end(), aper);
  
  // simultaneous fit on all objects
  SimFit zeFit(itSet->PhoRef);
  for_each(fidList.begin(), fidList.end(), simFit);

  cout << " LightCurveGuru::MonopolizeCPU() : releasing CPU \n";

}
