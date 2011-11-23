#ifndef IMAGESUBTRACTION__H
#define IMAGESUBTRACTION__H

#include "reducedimage.h"
#include "fileutils.h"
#include "kernelfitter.h"
#include "detection.h"

string SubtractedName(const string &RefName, const string &NewName);
    
//! For subtracting images using the Alard kernel fit technique
//! A basic assumption: Ref and New are already geometrically aligned
class ImageSubtraction : public ReducedImage, public KernelFitter {

  ReducedImageRef Ref, New;

 public :
  //! the constructor takes a copy of both ReducedImage.
  ImageSubtraction(const string &Name, const ReducedImageRef RefImage,
		   const ReducedImageRef NewImage, const bool NoSwap = false);
  
  //! as above but uses kernel found by other means. used for fakes SN
  ImageSubtraction(const string &Name, const ReducedImageRef RefImage,
		   const ReducedImageRef NewImage, const KernelFitter &Kfit);

  //! constuctor to read a produced subtraction
  ImageSubtraction(const string &Name);

  virtual const string  TypeName() const { return "ImageSubtraction"; }

  //! produces the subtracted image
  bool MakeFits() ;
  
  //! produces a weight with both ref and new counted for
  bool MakeWeight();
  
  //! carry out the detection. Default is matched filter plus aperture photometry.
  bool MakeCatalog();

  //! do cosmic but do not flag detection, only update weight
  bool MakeCosmic();

  // ! Dead frame is the OR of Ref and New dead frames.
  bool MakeDead();

  //!  Satur frame is the OR of Ref and New satur frames.
  bool MakeSatur();

  //!  Mask the saturation on the subtraction
  bool MaskSatur();

  //!  Mask null weights on the subtraction
  bool MaskNullWeight();
  
  string MatchedCatName() const { return Dir() + "matcheddet.list"; }
  string CatalogName() const { return Dir() + "det.list"; }
  string PsfCatalogName() const { return Dir() + "psfdet.list"; }
  string KernelName() const { return Dir() + "kernel.dat"; }

  //! useless clone
  ReducedImageRef Clone() const;
};

typedef CountedRef<ImageSubtraction> ImageSubtractionRef;

#endif // IMAGESUBTRACTION__H
