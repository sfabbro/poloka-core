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
  double refSeeing = Ref->GFSeeing();
  double newSeeing = New->GFSeeing();  
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

static string kernfitfile(const ReducedImageRef best, const ReducedImageRef worst) {
  ImagePair impair(best,worst);
  return KernelFitFile(impair) + ".dat";
}
void KernelFitter::WriteKernel(bool overwrite) const
{
  string kernfile = kernfitfile(worst,best);
  if (FileExists(kernfile)) {
    if (overwrite) remove(kernfile.c_str());
    else return;
  }
  cout << " Writing " << kernfile << endl;
  write(kernfile);
}

bool KernelFitter::ReadKernel()
{
  string kernfile = kernfitfile(worst,best);
  if (!FileExists(kernfile)) return false;
  cout << " Reading " << kernfile << endl;
  read(kernfile);
  return true;
}
