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
  bool writeMatrices; 
  bool subDirPerObject;
  
 public:
  
  LightCurveFile(const string &LCFileName);

  //! regular fit: fits objects in LCFile
  bool SimPhotFitAll(double vignette_size_n_seeing=-1) const;

  //! fit objects in calib catalog
  bool  SimPhotFitAllCalib(const string &CalibCatalog, const string &OutputCatalog, int itype=1, int Nmax=-1, double vignette_size_n_seeing=-1) const;

  const ReducedImage* GeomRef() const { return geomRef;}
  const ReducedImageList &Images() const { return images;}

  bool WriteVignettes() const { return writeVignettes;}
  bool WriteMatrices()  const { return writeMatrices;}
  void PleaseWriteVignettes() { writeVignettes = true;}
  void PleaseWriteMatrices() { writeMatrices = true;}

  bool OneDirPerObject() const { return subDirPerObject;}
  void PleaseOneDirPerObject() { subDirPerObject = true;}

  

};


#endif /* LIGHTCURVEFILE__H */
