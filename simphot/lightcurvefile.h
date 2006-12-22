#ifndef LIGHTCURVEFILE__H
#define LIGHTCURVEFILE__H


#include "objecttofit.h"
#include "reducedimage.h"


#include <list>
using namespace std;

class LightCurveFile
{
 private:
  ObjectToFitList objects;
  ReducedImageList images;
  ReducedImageRef geomRef;
  bool writeVignettes;
  bool subDirPerObject;
  string tupleFileName;
  
 public:
  
  LightCurveFile(const string &LCFileName);

  //! regular fit: fits objects in LCFile
  bool SimPhotFitAll() const;

  //! fit objects in calib catalog
  bool  SimPhotFitAllCalib(const string &CalibCatalog) const;

  const ReducedImage* GeomRef() const { return geomRef;}
  const ReducedImageList &Images() const { return images;}

  bool WriteVignettes() const { return writeVignettes;}
  void PleaseWriteVignettes() { writeVignettes = true;}

  bool OneDirPerObject() const { return subDirPerObject;}
  void PleaseOneDirPerObject() { subDirPerObject = true;}

  

};


#endif /* LIGHTCURVEFILE__H */
