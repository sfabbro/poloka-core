#ifndef TRANSFORMEDIMAGESADDFAKES__H
#define TRANSFORMEDIMAGESADDFAKES__H


#include "transformedimage.h"

class SENearStarList;
class Gtransfo;

class TransformedImageAddFakes : public TransformedImage {
 private:
  SENearStarList *NearStarRef;
  Gtransfo *One2Two;

 public:
  //Constructeurs
  TransformedImageAddFakes(const string &Name, 
			   const ReducedImage *Source, 
			   const ImageTransfo *transfo, 
			   SENearStarList *NSList, 
			   const Gtransfo *oneto2);
  
  bool MakeFits();
  ReducedImage *Clone() const;

  ~TransformedImageAddFakes();


};
#endif /* TRANSFORMEDIMAGESITHFAKES__H */
