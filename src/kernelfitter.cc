#include "kernelfitter.h"
#include "imagepair.h"
#include "polokaexception.h"






KernelFitter::KernelFitter(const ReducedImageRef &Ref, const ReducedImageRef &New, const bool NoSwap) 
{
  if (!(Ref->SamePhysicalSize(*New)))
    {
      cerr << " ERROR : Cannot subtract images of different sizes : " 
	   << Ref->Name() << " and " << New->Name() << endl;
      throw PolokaException(" Cannot fit kernel on images of different sizes ");
    }
  double refSeeing = Ref->Seeing();
  double newSeeing = New->Seeing();
  largestSeeing = (refSeeing > newSeeing) ? refSeeing : newSeeing;
  refIsBest =  (NoSwap || (refSeeing < newSeeing));
  if (refIsBest)
    {
      best = Ref;
      worst = New;
    }
  else
    {
      best = New;
      worst = Ref;
    }
  commonFrame = ImagePair(Ref,New).CommonFrame();
}


int KernelFitter::DoTheFit()
{
  ImagePair imPair(best,worst);
  return KernelFit::DoTheFit(imPair);  
  // clean up of images when deleting imPair;
}
