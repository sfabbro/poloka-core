// This may look like C code, but it is really -*- C++ -*-
#ifndef IMAGESUM__H
#define IMAGESUM__H


/* 

this is in principle not more than an image with the information needed
to sum it with other images 

*/

#include <vector>
#include <string>
#include "reducedimage.h"

enum StackingMethod { WeightedAverage = 1, 
		      ClippedWeightedAverage = 2, 
		      Median = 3,
		      SUnSet = 4};

enum WeightingMethod { PointSourceOptimal = 1, 
		       ExtendedSourceOptimal = 2, 
		       NoGlobalWeighting = 3, //Local weights used, global=1
		       NoWeightsAtAll = 4, //Local and global weights ignored
		       WUnSet = 5};

string name_of_stackingMethod(const StackingMethod stackingMethod);
string name_of_weightingMethod(const WeightingMethod weightingMethod);

class FitsHeader;


//! store the necessary components for weighting and computing the stacking
class Component {
public :
  double backVar; // sigma^2(back)
  double seeing;
  double photomRatio; // (this image) / (photometric reference)
  double globalWeight;
  double averageWeight; // used only to weight the various header words.
  bool hasWeight, hasDead, hasSatur, hasCosmic;
  string reducedImageName;
  ReducedImageRef Ri; //! not to be stored
  
  Component(ReducedImage *RI, const double &PhotomRatio, 
	    const WeightingMethod weightingMethod);
  
  
  //typedef vector<Component>::iterator ComponentIterator;
  //typedef vector<Component>::const_iterator ComponentCIterator;
  void dump(ostream& s = cout) const;

  ~Component();
private:
  void SetGlobalWeights(const WeightingMethod weightingMethod);
};



//! A class that handles coadding. Shift and coadd is handled by 
//ImagesAlignAndSum()
/*! Shift is performed using the TransformedImage class. Addition is performed
  using ImageSum. Depending on the numner of images involved, we perform either
  an actual sum, or a "clipped mean".
*/
class ImageSum : public ReducedImage {
   void init();
private:

  vector<Component> components;
  StackingMethod stackingMethod;
  WeightingMethod weightingMethod;
 
 double zero_point_ref;

  double totalWeight;
  double seeing;
  double backLevel;
  double sigmaBack;
  double saturation;
  double ActualNoise;
  double originalsatur;
  double exposure;
  double readnoise;
  double flatnoise;
  double julianDate;
  double airmass;
  double originalskylevel;
  string photomReferenceName;

  Frame intersection;
  
  
public :
  ImageSum(const string &Name, ReducedImageList &Images,
	   const ReducedImage *PhotomReference = NULL,
	   const WeightingMethod AWMethod = WUnSet, 
	   const StackingMethod ASMethod = SUnSet);

  ImageSum(const string &Name);
  ImageSum(){};
  ReducedImageList Components() const;

  virtual const string  TypeName() const { return "ImageSum";}

  bool MakeFits() ;
  bool MakeCatalog() ; 
  bool MakeDead();
  bool MakeSatur();
  bool MakeWeight();
  ReducedImage *Clone() const;

  void dump(ostream & s = cout) const;

  typedef vector<Component>::iterator ComponentIterator;
  typedef vector<Component>::const_iterator ComponentCIterator;

  ~ImageSum();

private:
  void FitsHeaderFill();

};


//! aligns and sum in a single call.

ImageSum* ImagesAlignAndSum(const ReducedImageList &ToSum, 
			    const ReducedImage &GeomReference, 
			    const string &SumName, const int ToDo,
			    const ReducedImage *PhotomRef = NULL,
			    const WeightingMethod AWMethod = WUnSet, 
			    const StackingMethod ASMethod = SUnSet);

//!
bool ImagesAlignAndSum(const vector<string> &ToSum, 
		       const string &Reference, 
		       const string &SumName, const int ToDo);

#endif /* IMAGESUM__H */
