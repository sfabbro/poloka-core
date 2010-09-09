#ifndef IMAGEPAIR__H
#define IMAGEPAIR__H

#include <string>
#include "frame.h"
#include "countedref.h"

class ReducedImage;
class FitsImage;
class Image;

class ImagePair {

  typedef CountedRef<ReducedImage> ReducedImageRef;

  ReducedImageRef best,worst;
  double bestSeeing,worstSeeing;
  double bestGain, worstGain;
  Frame commonFrame;
  const Image* images[4];

 public :
  ImagePair(const ReducedImageRef &Best, const ReducedImageRef &Worst);

  const ReducedImage* Best() const { return best;}
  const ReducedImage* Worst() const { return worst;}

  //! Which = "Best","Worst", BestWeight", WorstWeight"
  const Image&  BestImage();
  const Image&  BestWeight();
  const Image&  WorstImage();
  const Image&  WorstWeight();
  
  double BestSeeing() const { return bestSeeing;}
  double WorstSeeing() const {return worstSeeing;}
  double BestGain() const {return bestGain;}
  double WorstGain() const {return worstGain;}

  const Frame& CommonFrame() const { return commonFrame;}

  ~ImagePair();

 private :
// forbid copies : pointers on images should not be deleted twice
  ImagePair(const ImagePair &); 
  void operator = (const ImagePair &Right); 

};


#endif /*IMAGEPAIR__H */
